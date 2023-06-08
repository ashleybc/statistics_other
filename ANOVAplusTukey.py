import pandas as pd
from scipy import stats
import numpy as np
from statsmodels.stats.multicomp import pairwise_tukeyhsd

root="/Users/ashley/Dropbox/Research/Graduate/SoMAS/GordonsLab/FGL/code/"
mydf=pd.read_csv(root+"inputs/GL5_enr_data.csv")
mydf=mydf[~mydf["treatment"].str.contains("Poly")]
mydf.reset_index(inplace=True)

#just constructing a full name id out of depth and treatment in case i want to use later
newNames=[]
for count,t in enumerate(mydf["treatment"]):
    #print(count,t,mydf["Depth"][count])
    full_name=t+"_"+str(mydf["Depth"][count])
    print(full_name)
    newNames.append(full_name)
    
mydf["Full_ID"]=newNames

def one_way_ANOVA_pre_proc(df,group_col,val_col):

    #setting up the data structure for one way ANOVA by creating a dictionary from treatment-measurement pairs
    key_val_pairs=[]

    #using pandas groupby function returns tuple as (group,data)
    for g in df.groupby(group_col):
        #print(g)
        #first item is treatment (group)
        group=g[0]
        data=list(g[1][val_col].values)
        #second item is the specified data. use .values to get out of dataframe into array or list format
        #print(data)
        #create key,val tuple and append to a list
        key_valpair=(group,data)
        #print(key_valpair)
        key_val_pairs.append(key_valpair)

    #create the dictionary from the list of tuples
    keyval_Dict=dict(key_val_pairs)
    
    return keyval_Dict

#testing one-way ANOVA on treatment only (disregarding depth)
test_dict=one_way_ANOVA_pre_proc(mydf,"treatment","IC assimilation rate (umole IC/L/d)")

#now in the one way ANOVA, can call any groupings by dictionary key
# stats f_oneway functions takes the groups as input and returns F and P-value
fvalue, pvalue = stats.f_oneway(test_dict['C14+S2O3+NO3'],
                                test_dict['C14+S2O3'],
                                test_dict['Dark'])
print(fvalue,pvalue)
# perform multiple pairwise comparison (Tukey HSD)
m_comp = pairwise_tukeyhsd(endog=mydf['volume'], groups=mydf['engine'], alpha=0.05)
print(m_comp)

#unique depths
unique_dps=list(mydf.Depth.unique())


##setting up the data structure for one way ANOVA by creating a dictionary from treatment-measurement pairs
#key_val_pairs=[]
#
##using pandas groupby function returns tuple as (group,data)
#for g in mydf.groupby("treatment"):
#    print(g)
#    #first item is treatment (group)
#    group=g[0]
#    data=list(g[1]["IC assimilation rate (umole IC/L/d)"].values)
#    #second item is the specified data. use .values to get out of dataframe into array or list format
#    print(data)
#    #create key,val tuple and append to a list
#    key_valpair=(group,data)
#    print(key_valpair)
#    key_val_pairs.append(key_valpair)
#
##create the dictionary from the list of tuples
#myDict=dict(key_val_pairs)


#if want to specify treatments
#for g in mydf.groupby("Depth"):
#    depth=g[0]
#    print("Depth:",depth)
#    print(g[1],"\n")
#    #testing one-way ANOVA on treatment only (disregarding depth)
#    test_dict=one_way_ANOVA_pre_proc(g[1],"treatment","IC assimilation rate (umole IC/L/d)")
#
#    print(test_dict,"\n")
#    #now in the one way ANOVA, can call any groupings by dictionary key
#    # stats f_oneway functions takes the groups as input and returns F and P-value
#    fvalue, pvalue = stats.f_oneway(test_dict['C14+S2O3'],test_dict['Dark'])
#    print(fvalue,pvalue,"\n\n")


#do available treatments per depth:
Depth=[]
pvals=[]
fvals=[]
Treatments=[]

for g in mydf.groupby("Depth"):
    depth=g[0]
    print("Depth:",depth)
    Depth.append(depth)
    print(g[1],"\n")
    #pre-processing for one way ANOVA. need a dictionary of treatments and corresponding response values
    test_dict=one_way_ANOVA_pre_proc(g[1],"treatment","IC assimilation rate (umole IC/L/d)")

    print(test_dict,"\n")
    avail_keys=list(test_dict.keys())
    print("available keys:",avail_keys)

    #Because not all keys (treatments) were done per depth, need to determine available keys and vals
    keyList=[]
    valList=[]

    for k in avail_keys:
        keyList.append(k)
        valList.append(test_dict[k])

    print(valList)

    Treatments.append(";".join(keyList))
    #now in the one way ANOVA, can call any groupings by dictionary key
    # stats f_oneway functions takes the groups as input and returns F and P-value
    #use * to unpack nested list of lists
    fvalue, pvalue = stats.f_oneway(*valList)
    pvals.append(pvalue)
    fvals.append(fvalue)
    print("f,p",fvalue,pvalue,keyList)

    #screen, because if only one treatment, cannot perform pairwise comparison
    if len(keyList)>1:
        #response variable is ICA, grouping in this case is treatment
        m_comp=pairwise_tukeyhsd(endog=g[1]["IC assimilation rate (umole IC/L/d)"],groups=g[1]["treatment"],alpha=0.05)
        print("tukey:",m_comp,"\n\n")
        tukey_df = pd.DataFrame(data=m_comp._results_table.data[1:], columns=m_comp._results_table.data[0])
        tukey_df.name=str(depth)+"_GL5enr_tukey"
        tukey_df.to_csv(root+"outputs/"+tukey_df.name+".csv",index=False)

#once iterated over all depths, make one ANOVA df and export
anovaDF=pd.DataFrame({"Depth (m)":Depth,"f":fvals,"p":pvals,"Treatment":Treatments})
anovaDF.set_index("Depth (m)",inplace=True)
anovaDF.to_csv(root+"outputs/GL5enrichment_incs_onewayanovatests.csv")
