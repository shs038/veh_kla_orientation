```
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
from collections import Counter
%matplotlib inline
```
#read in data
```
motif_orientation_df = pd.read_csv('/home/jtao/for_shengnan/motif_strand_C57BL6J.tsv', sep='\t')
motif_orientation_df.index = motif_orientation_df['ID'].values
del motif_orientation_df['ID']
```
#vehicle orientation
```
veh_df=motif_orientation_df[motif_orientation_df['Factors'].str.contains('c57bl6_atac_veh')]
del veh_df['Factors']
del veh_df['chr']
```
#KLA orientation 
```
kla_df=motif_orientation_df[motif_orientation_df['Factors'].str.contains('c57bl6_atac_kla')]
del kla_df['Factors']
del kla_df['chr']
```
#function for counting orientations
```
def count_orientation(motif_orientation):
    *input:a pandas dataframe contains motifs orientation data*
    *output:a 3 rows pandas dataframe contains counts of each orientation of each motifs*
    motifs = motif_orientation.columns.values #save motifs identity 
    *create a zeros matrix to sotre future orientation count data*
    zero_data = np.zeros((3,motif_orientation.shape[1]),dtype=np.int)
    *convert zeros matrix to zeros dataframe*
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=['+','-','?']
    *loop to count the orientations of all motifs*
    for i in range (motif_orientation.shape[1]):
        one_motif=motif_orientation.ix[:,i].values*retrieve orientations of one motif*
        one_motif=list(one_motif)*convert to list*
        c=Counter(one_motif)*count each orientation*
        c=dict(c)*convert to dictionary*
        c=pd.DataFrame.from_dict(c,orient='index')*convert to dataframe*
        count_frame.ix[:,i]=c.ix[:,0]*store orientation count in zeros dataframe*
    return count_frame
```
#count orientation
```
veh_orientation=count_orientation(veh_df)
kla_orientation=count_orientation(kla_df)
```
#transpose the count dataframe 
```
veh_orientation = veh_orientation.T
kla_orientation = kla_orientation.T
```
#count how many times two motifs that co-occur with each other both have sense orientation 
```
def count_both_sense(motif_orientation):
    *input:a pandas dataframe contains motifs orientation data*
    *output:a pandas dataframe contains how many times each pair of motifs both have + orientation.*
    motifs = motif_orientation.columns.values *save motifs identity*
    *creaty a zeros matrix to sotre future orientation count data*
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    *conver zeros matrix to zeros dataframe*
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        *find the loci where the motif occur with sense orientation*
        logical_col_i=motif_orientation.ix[:,i]=='+' 
        for j in range (i+1,motif_orientation.shape[1]):
            #find the loci where the motif occur with sense orientation
            logical_col_j=motif_orientation.ix[:,j]=='+'
            #find the loci where both of the motifs occur with sense orientation
            logical_input=1*logical_col_i+1*logical_col_j==2
            #count how many times both of the motifs occur with sense orientation
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    #reshape dataframe
    Pairs=[]
    Count=[]
    #loop in part of count data that contain meaning counting
    for i in range (count_frame.shape[1]-1):
        for j in range (i+1,count_frame.shape[1]):
            #put motif pair and correlation into the empty list
            motif_pairs=(motifs[i],motifs[j])
            Pairs.append(motif_pairs)
            Count.append(count_frame.ix[i,j])
    #reshape the dataframe
    reshaped_frame = pd.DataFrame({'Count': Count}, index=Pairs)
    reshaped_frame['Orientation']='+/+' #add orientation to the dataframe
    reshaped_frame=reshaped_frame[['Orientation','Count']]# put Orientation in front of Count
    return reshaped_frame
```

#count how many times two motifs that co-occur with each other both have antisense orientation 
```
def count_both_antisense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs both have 
            - orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #create a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        logical_col_i=motif_orientation.ix[:,i]=='-'
        for j in range (i+1,motif_orientation.shape[1]):
            logical_col_j=motif_orientation.ix[:,j]=='-'
            logical_input=1*logical_col_i+1*logical_col_j==2
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    #reshape dataframe
    Pairs=[]
    Count=[]
    #loop in part of count data that contain meaning counting
    for i in range (count_frame.shape[1]-1):
        for j in range (i+1,count_frame.shape[1]):
            #put motif pair and correlation into the empty list
            motif_pairs=(motifs[i],motifs[j])
            Pairs.append(motif_pairs)
            Count.append(count_frame.ix[i,j])
    #reshape the dataframe
    reshaped_frame = pd.DataFrame({'Count': Count}, index=Pairs)
    reshaped_frame['Orientation']='-/-' #add orientation to the dataframe
    reshaped_frame=reshaped_frame[['Orientation','Count']]# put Orientation in front of Count
    return reshaped_frame
```
#count how many times two motifs co-occur with each other have sense/antisense orientation 
```
def count_sense_antisense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs  have 
            +/- orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #create a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        logical_col_i=motif_orientation.ix[:,i]=='+'
        for j in range (i+1,motif_orientation.shape[1]):
            logical_col_j=motif_orientation.ix[:,j]=='-'
            logical_input=1*logical_col_i+1*logical_col_j==2
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    #reshape dataframe
    Pairs=[]
    Count=[]
    #loop in part of count data that contain meaning counting
    for i in range (count_frame.shape[1]-1):
        for j in range (i+1,count_frame.shape[1]):
            #put motif pair and correlation into the empty list
            motif_pairs=(motifs[i],motifs[j])
            Pairs.append(motif_pairs)
            Count.append(count_frame.ix[i,j])
    #reshape the dataframe
    reshaped_frame = pd.DataFrame({'Count': Count}, index=Pairs)
    reshaped_frame['Orientation']='+/-' #add orientation to the dataframe
    reshaped_frame=reshaped_frame[['Orientation','Count']]# put Orientation in front of Count
    return reshaped_frame
```
#count how many times two motifs co-occur with each other have antisense/sense orientation 
```
def count_antisense_sense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs  have 
            -/+ orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #creaty a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        logical_col_i=motif_orientation.ix[:,i]=='-'
        for j in range (i+1,motif_orientation.shape[1]):
            logical_col_j=motif_orientation.ix[:,j]=='+'
            logical_input=1*logical_col_i+1*logical_col_j==2
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    #reshape dataframe
    Pairs=[]
    Count=[]
    #loop in part of count data that contain meaning counting
    for i in range (count_frame.shape[1]-1):
        for j in range (i+1,count_frame.shape[1]):
            #put motif pair and correlation into the empty list
            motif_pairs=(motifs[i],motifs[j])
            Pairs.append(motif_pairs)
            Count.append(count_frame.ix[i,j])
    #reshape the dataframe
    reshaped_frame = pd.DataFrame({'Count': Count}, index=Pairs)
    reshaped_frame['Orientation']='-/+' #add orientation to the dataframe
    reshaped_frame=reshaped_frame[['Orientation','Count']]# put Orientation in front of Count
    return reshaped_frame
```
#dataframe of veh orientation count
```
veh_sense_sense=count_both_sense(veh_df)
veh_antisense_antisense=count_both_antisense(veh_df)
veh_sense_antisense=count_sense_antisense(veh_df)
veh_antisense_sense=count_antisense_sense(veh_df)
```
#concatenate data for all four orientation
```
veh_frames = [veh_sense_sense,veh_antisense_antisense,
          veh_sense_antisense,veh_antisense_sense]
veh_cooccur_orientation = pd.concat(veh_frames)
```
#dataframe of kla orientation count
```
kla_sense_sense=count_both_sense(kla_df)
kla_antisense_antisense=count_both_antisense(kla_df)
kla_sense_antisense=count_sense_antisense(kla_df)
kla_antisense_sense=count_antisense_sense(kla_df)
```
#concatenate data for all four orientation
```
kla_frames = [kla_sense_sense,kla_antisense_antisense,
          kla_sense_antisense,kla_antisense_sense]
kla_cooccur_orientation = pd.concat(kla_frames)
```
#normalize orientation count for veh
```
veh_df_add=veh_sense_sense['Count'].add(veh_antisense_antisense['Count'], fill_value=0)
veh_df_add=veh_df_add.add(veh_sense_antisense['Count'], fill_value=0)
veh_df_add=veh_df_add.add(veh_antisense_sense['Count'], fill_value=0)
veh_cooccur_orientation['Total']=veh_df_add
veh_cooccur_orientation = veh_cooccur_orientation[veh_cooccur_orientation.Total != 0]
veh_cooccur_orientation['Normalized_Count']=veh_cooccur_orientation['Count']/veh_cooccur_orientation['Total']
```
#normalize orientation count for kla
```
kla_df_add=kla_sense_sense['Count'].add(kla_antisense_antisense['Count'], fill_value=0)
kla_df_add=kla_df_add.add(kla_sense_antisense['Count'], fill_value=0)
kla_df_add=kla_df_add.add(kla_antisense_sense['Count'], fill_value=0)
kla_cooccur_orientation['Total']=kla_df_add
kla_cooccur_orientation = kla_cooccur_orientation[kla_cooccur_orientation.Total != 0]
kla_cooccur_orientation['Normalized_Count']=kla_cooccur_orientation['Count']/kla_cooccur_orientation['Total']
```
#subset of veh_cooccur_orientation for each orientation pair without Normalized_Count=0
```
veh_cooccur_sense_sense=veh_cooccur_orientation[veh_cooccur_orientation['Orientation']=='+/+']
veh_cooccur_antisense_antisense=veh_cooccur_orientation[veh_cooccur_orientation['Orientation']=='-/-']
veh_cooccur_sense_antisense=veh_cooccur_orientation[veh_cooccur_orientation['Orientation']=='+/-']
veh_cooccur_antisense_sense=veh_cooccur_orientation[veh_cooccur_orientation['Orientation']=='-/+']
```
#subset of kla_cooccur_orientation for each orientation pair without Normalized_Count=0
```
kla_cooccur_sense_sense=kla_cooccur_orientation[kla_cooccur_orientation['Orientation']=='+/+']
kla_cooccur_antisense_antisense=kla_cooccur_orientation[kla_cooccur_orientation['Orientation']=='-/-']
kla_cooccur_sense_antisense=kla_cooccur_orientation[kla_cooccur_orientation['Orientation']=='+/-']
kla_cooccur_antisense_sense=kla_cooccur_orientation[kla_cooccur_orientation['Orientation']=='-/+']
```
#veh orientation array for chi square test
```
motifpair=veh_cooccur_sense_sense.index.values
vss=veh_cooccur_sense_sense.ix[:,1].values
vsa=veh_cooccur_sense_antisense.ix[:,1].values
vas=veh_cooccur_antisense_sense.ix[:,1].values
vaa=veh_cooccur_antisense_antisense.ix[:,1].values
veh_all_orientation=np.array([vss,vsa,vas,vaa])
veh_all_orientation=veh_all_orientation.T
```
#kla orientation array for chi square test
```
kss=kla_cooccur_sense_sense.ix[:,1].values
ksa=kla_cooccur_sense_antisense.ix[:,1].values
kas=kla_cooccur_antisense_sense.ix[:,1].values
kaa=kla_cooccur_antisense_antisense.ix[:,1].values
kla_all_orientation=np.array([kss,ksa,kas,kaa])
kla_all_orientation=kla_all_orientation.T
```
# Given a set of loci L where motif J is in a given orientation O, does that subset of L where 
## motf Z is present ahve a bias towards +/-
## Fisher's  exact test 
## Null hypothesis: orientation of motif J and that of motif Z are independent at loci L.
```
def Chisquare_test(orientation):
    Chi_array=np.array(np.zeros((1,4),dtype=np.int))
    #store
    P_value=[]
    for i in range(len(orientation)):
        Chi_array=orientation[i,:]
        chisq,p=stats.chisquare(Chi_array)
        P_value.append(p)
    P_df=pd.DataFrame({'P_value': P_value}, index=motifpair)
    return P_df
```
# Chi squiare test of orientation, corrected by times of compare
```
veh_p=Chisquare_test(veh_all_orientation)*195
kla_p=Chisquare_test(kla_all_orientation)*195
```
# dataframe for p_value of motif pairs orientation 
```
orientation_p_df=pd.concat([veh_p, kla_p], axis=1)
orientation_p_df.columns = ['veh', 'kla']
orientation_p_df['veh-kla']=orientation_p_df['veh']-orientation_p_df['kla']
orientation_p_df.sort('veh-kla',ascending=False)
```
#plot difference
```
sns.distplot(orientation_p_df['vel-kla'])
plt.ylabel('Frequency')
plt.xlabel('vel-kla')
plt.title('orientation p_difference under vel and kal') 
```
