
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

not_plot_upper_right_line=1;


def individual_clone_dynamics_sampling_M1(P_array,sample_1):
    # random sample from the same set in both T1 and T2

    #P_ex=0.7; # killing mouse, get HSC via FACS
    #P_sp: first partitioning: one for inDrops; the other for secondary transplantation
    #P_nd=0.7; # Non Dropout for the first inDrops
    #p2=0.49; # HSC extraction + Non Dropout rate for the second inDrops
    
    
    #alpha_1=10; #clonal expansion rate for the first transplantation
    #alpha_2=1;  #clonal expansion rate for the secondary transplantation
    #beta=0.5; # chance to be successfully engrafted
    P_ex=P_array[0];
    P_indrops=P_array[1];
    P_sp=P_array[2];
    P_nd=P_array[3];
    beta=P_array[4];


    # HSC expansion during the first transplantation
    m1_0=np.random.choice(sample_1);

    # kill the mouse, extact the bone marrow cells plus FACS
    m_star_0=np.random.binomial(m1_0,P_ex,1)[0];
    
    # split the population into two parts, one for inDrops (m_seq), and the other for secondary transplantation (m_star)
    m_seq=np.random.binomial(m_star_0,P_indrops,1)[0]; # the fraction is fixed to be 50%

    # first indrops result
    m1=np.random.binomial(m_seq,P_nd,1)[0];

    # the other part
    m_star=np.random.binomial(m_star_0-m_seq,P_sp,1)[0];  # normally P_sp=1, if divide agagin, P_sp=0.5;

    #engraftment of the HSC during secondary transplantation
    m_en=np.random.binomial(m_star,beta,1)[0];

    #clonal expansion for the second transplantation
    m2_0=np.sum([np.random.choice(sample_1) for j in range(m_en)]);

    # observed clone size during the second inDrops
    m2=np.random.binomial(m2_0,P_ex*P_nd,1)[0];
    return [int(m1),int(m2)]
    


def individual_clone_dynamics_sampling_M1_Kit1(P_array,HSC_list,Kit_list):
    # random sample from the same set (include both HSC and Kit information) in both T1 and T2

    # for each clonal expansion, output now only HSC clone size but also Kit+ clone size
    
    # HSC expansion during the first transplantation

    P_ex=P_array[0];
    P_indrops=P_array[1];
    P_sp=P_array[2];
    P_nd=P_array[3];
    beta=P_array[4];
    
    index_list=np.array(range(len(HSC_list)),dtype="int");
    index_1=np.random.choice(index_list);
    m1_0=HSC_list[index_1];
    Kit1_0=Kit_list[index_1];

    # kill the mouse, extact the bone marrow cells plus FACS
    m_star_0=np.random.binomial(m1_0,P_ex,1)[0];
    Kit_star_0=np.random.binomial(Kit1_0,P_ex,1)[0];
    
    # split the population into two parts, one for inDrops (m_seq), and the other for secondary transplantation (m_star)

    m_seq=np.random.binomial(m_star_0,P_indrops,1)[0]; # the fraction is fixed to be 50%
    Kit_seq=np.random.binomial(Kit_star_0,P_indrops,1)[0]; 

    # first indrops result
    m1=np.random.binomial(m_seq,P_nd,1)[0];
    Kit1=np.random.binomial(Kit_seq,P_nd,1)[0];


    m_star=np.random.binomial(m_star_0-m_seq,P_sp,1)[0]; # normally P_sp=1, if divide agagin, P_sp=0.5;
#    print("m_star_0={}, m_seq={}, m_star={}".format(m_star_0,m_seq,m_star))

    #engraftment of the HSC during secondary transplantation
    m_en=np.random.binomial(m_star,beta,1)[0];

    #clonal expansion for the second transplantation
    index_2=[np.random.choice(index_list) for j in range(m_en)];
    #index_2=index_2.astype(int);
    m2_0=np.sum(HSC_list[index_2]);
    Kit2_0=np.sum(Kit_list[index_2]);

    # observed clone size during the second inDrops
    m2=np.random.binomial(m2_0,P_ex*P_nd,1)[0];
    Kit2=np.random.binomial(Kit2_0,P_ex*P_nd,1)[0];
    return [int(m1),int(m2),int(Kit1),int(Kit2)]





def individual_clone_dynamics_sampling_M1_Kit1_two_mouse(P_array,HSC_list,Kit_list):
    # random sample from the same set (include both HSC and Kit information) in both T1 and T2

    # for each clonal expansion, output now only HSC clone size but also Kit+ clone size
    
    # HSC expansion during the first transplantation

    P_ex=P_array[0];
    P_indrops=P_array[1];
    P_sp=P_array[2];
    P_nd=P_array[3];
    beta=P_array[4];
    

    index_list=np.array(range(len(HSC_list)),dtype="int");
    index_1=np.random.choice(index_list);
    m1_0=HSC_list[index_1];
    Kit1_0=Kit_list[index_1];

    # kill the mouse, extact the bone marrow cells plus FACS
    m_star_0=np.random.binomial(m1_0,P_ex,1)[0];
    Kit_star_0=np.random.binomial(Kit1_0,P_ex,1)[0];
    
    # split the population into two parts, one for inDrops (m_seq), and the other for secondary transplantation (m_star)

    m_seq=np.random.binomial(m_star_0,P_indrops,1)[0]; # the fraction is fixed to be 50%
    Kit_seq=np.random.binomial(Kit_star_0,P_indrops,1)[0]; 

    # first indrops result
    m1=np.random.binomial(m_seq,P_nd,1)[0];
    Kit1=np.random.binomial(Kit_seq,P_nd,1)[0];

    m_star_mouse_1=np.random.binomial(m_star_0-m_seq,P_sp,1)[0]; # normally P_sp=1, if divide agagin, P_sp=0.5;
    m_star_mouse_2=m_star_0-m_seq-m_star_mouse_1    

    #engraftment of the HSC during secondary transplantation
    m_en_mouse_1=np.random.binomial(m_star_mouse_1,beta,1)[0];
    m_en_mouse_2=np.random.binomial(m_star_mouse_2,beta,1)[0];

    #clonal expansion for the first mouse in 2T 
    index_2_mouse_1=[np.random.choice(index_list) for j in range(m_en_mouse_1)];
    m2_0_mouse_1=np.sum(HSC_list[index_2_mouse_1]);
    Kit2_0_mouse_1=np.sum(Kit_list[index_2_mouse_1]);

    #clonal expansion for the second mouse in 2T 
    index_2_mouse_2=[np.random.choice(index_list) for j in range(m_en_mouse_2)];
    m2_0_mouse_2=np.sum(HSC_list[index_2_mouse_2]);
    Kit2_0_mouse_2=np.sum(Kit_list[index_2_mouse_2]);

    # observed clone size during the second inDrops
    m2_mouse_1=np.random.binomial(m2_0_mouse_1,P_ex*P_nd,1)[0];
    m2_mouse_2=np.random.binomial(m2_0_mouse_2,P_ex*P_nd,1)[0];
    Kit2_mouse_1=np.random.binomial(Kit2_0_mouse_1,P_ex*P_nd,1)[0];
    Kit2_mouse_2=np.random.binomial(Kit2_0_mouse_2,P_ex*P_nd,1)[0];
    return [int(m1),int(m2_mouse_1),int(m2_mouse_2),int(Kit1),int(Kit2_mouse_1),int(Kit2_mouse_2)]





def individual_clone_dynamics_sampling_M1_Kit1_switching(P_array,expan_thres,activ_thres,HSC_list,Kit_list):
    # random sample from the same set (include both HSC and Kit information) in both T1 and T2

    # for each clonal expansion, output now only HSC clone size but also Kit+ clone size
    
    # HSC expansion during the first transplantation

    P_ex=P_array[0];
    P_indrops=P_array[1];
    P_sp=P_array[2];
    P_nd=P_array[3];
    beta=P_array[4];
    
    low_expansion_list_HSC=HSC_list[HSC_list<=expan_thres];
    low_expansion_list_Kit=Kit_list[HSC_list<=expan_thres];
    
    index_list=np.array(range(len(HSC_list)),dtype="int");
    index_1=np.random.choice(index_list);
    m1_0=HSC_list[index_1];
    Kit1_0=Kit_list[index_1];
    
    # define the sample list based on the activity in 1T
    # default sample list
    second_sample_list_HSC=HSC_list;
    second_sample_list_Kit=Kit_list;
    second_index_list=index_list;

        
#             #if m1_0>1/(P_ex*P_sp*P_nd): #this step, filter the noise.  Only when sufficiently large clone, the activity is accurate
#     #print("m1_0",m1_0);
    activity=-1;
    if m1_0>1:
        pseudo=1; # seems to be unimportant
        activity=(Kit1_0+pseudo)/(m1_0+pseudo);
        if activity>activ_thres:
            #print("activity",activity)
            second_sample_list_HSC=low_expansion_list_HSC;
            second_sample_list_Kit=low_expansion_list_Kit;
            second_index_list=np.array(range(len(second_sample_list_HSC)),dtype="int");
            

    # kill the mouse, extact the bone marrow cells plus FACS
    m_star_0=np.random.binomial(m1_0,P_ex,1)[0];
    Kit_star_0=np.random.binomial(Kit1_0,P_ex,1)[0];
    
    # split the population into two parts, one for inDrops (m_seq), and the other for secondary transplantation (m_star)
    m_seq=np.random.binomial(m_star_0,P_indrops,1)[0]; # the fraction is fixed to be 50%
    Kit_seq=np.random.binomial(Kit_star_0,P_indrops,1)[0]; 

    # first indrops result
    m1=np.random.binomial(m_seq,P_nd,1)[0];
    Kit1=np.random.binomial(Kit_seq,P_nd,1)[0];

    m_star=np.random.binomial(m_star_0-m_seq,P_sp,1)[0]; # normally P_sp=1, if divide agagin, P_sp=0.5;
    
    
    #if m1_0>1/(P_ex*P_sp*P_nd): #this step, filter the noise.  Only when sufficiently large clone, the activity is accurate
    #print("m1_0",m1_0);
#     activity=-1;
#     if m1>0:
#         pseudo=1; # seems to be unimportant
#         activity=(Kit1+pseudo)/(m1+pseudo);
#         if activity>activ_thres:
#             #print("activity",activity)
#             second_sample_list_HSC=low_expansion_list_HSC;
#             second_sample_list_Kit=low_expansion_list_Kit;
#             second_index_list=np.array(range(len(second_sample_list_HSC)),dtype="int");

    #engraftment of the HSC during secondary transplantation
    m_en=np.random.binomial(m_star,beta,1)[0];

    #clonal expansion for the second transplantation
    index_2=[np.random.choice(second_index_list) for j in range(m_en)];
    #index_2=index_2.astype(int);
    m2_0=np.sum(second_sample_list_HSC[index_2]);
    Kit2_0=np.sum(second_sample_list_HSC[index_2]);


    # observed clone size during the second inDrops
    m2=np.random.binomial(m2_0,P_ex*P_nd,1)[0];
    Kit2=np.random.binomial(Kit2_0,P_ex*P_nd,1)[0];
    
        
    return [int(m1),int(m2),int(Kit1),int(Kit2)]


def multi_clone_model_sampling_M1(N_initial_clone,P_array,inferred_M1):
    '''The parameter-free model, only ONE kind of HSC, each cell generate m progenies by randomly sampling from M1¶'''

    P_ex=P_array[0];
    P_sp=P_array[1];
    P_nd=P_array[2];
    beta=P_array[3];

    N=N_initial_clone;
    data=np.zeros((N,3));
#     print("Total number of initial clone:",N)
#     print("Enfrating rate:",beta);
#     print("Stem cell extraction rate:",P_ex)
#     print("Nondrop out rate:", P_nd)

    for j in range(0,N):
        [m1,m2]=individual_clone_dynamics_sampling_M1(beta,P_ex,P_sp,P_nd,inferred_M1);
        data[j,0]=m1;
        data[j,1]=m2;
        if m2==0 and m1==0:
            data[j,2]=0;
        elif m2>0 and m1==0:
            data[j,2]=0;
        else:
            data[j,2]=(m2)/(m1);
            
    return data;


def multi_clone_model_sampling_M1_kit(N_initial_clone,P_array,inferred_M1,inferred_Kit1):
    '''another parameter-free model, only ONE kind of HSC, with also Kit positive cells'''


    N=N_initial_clone;
    data=np.zeros((N,4));
#     print("Total number of initial clone:",N)
#     print("Enfrating rate:",beta);
#     print("Stem cell extraction rate:",P_ex)
#     print("Nondrop out rate:", P_nd)

    for j in range(0,N):
        [m1,m2,Kit1,Kit2]=individual_clone_dynamics_sampling_M1_Kit1(P_array,inferred_M1,inferred_Kit1);
        data[j,0]=m1;
        data[j,1]=m2;
        data[j,2]=Kit1;
        data[j,3]=Kit2;
            
    return data;



def multi_clone_model_sampling_M1_kit_two_mouse(N_initial_clone,P_array,inferred_M1,inferred_Kit1):
    '''another parameter-free model, only ONE kind of HSC, with also Kit positive cells. Here, we assume the other half of HSCs are further divided into half, each going to a different mouse. We are interested in the amount of correlation generated by the null model'''


    N=N_initial_clone;
    data=np.zeros((N,6));
#     print("Total number of initial clone:",N)
#     print("Enfrating rate:",beta);
#     print("Stem cell extraction rate:",P_ex)
#     print("Nondrop out rate:", P_nd)

    for j in range(0,N):
        [m1,m2_1,m2_2,Kit1,Kit2_1,Kit2_2]=individual_clone_dynamics_sampling_M1_Kit1_two_mouse(P_array,inferred_M1,inferred_Kit1);
        data[j,0]=m1;
        data[j,1]=m2_1;
        data[j,2]=m2_2;
        data[j,3]=Kit1;
        data[j,4]=Kit2_1;
        data[j,5]=Kit2_2;
            
    return data;




def multi_clone_model_activity_based_sampling_M1_Kit(N_initial_clone,P_array,expan_thres,activ_thres,inferred_M1,inferred_Kit1):

    
    N=N_initial_clone;
    #low_threshold=5; # threshold to define low-expansion list
    data=np.zeros((N,4));
#     print("Total number of initial clone:",N)
#     print("Enfrating rate:",beta);
#     print("Stem cell extraction rate:",P_ex)
#     print("Nondrop out rate:", P_nd)

    for j in range(0,N):
        [m1,m2,Kit1,Kit2]=individual_clone_dynamics_sampling_M1_Kit1_switching(P_array,expan_thres,activ_thres,inferred_M1,inferred_Kit1);
        data[j,0]=m1;
        data[j,1]=m2;
        data[j,2]=Kit1;
        data[j,3]=Kit2;
            
    return data;



def plot_average_expansion(x,y,my_list,plot_fig,fig_dir):
    ''' compute the local average of y according to the predefined clustering based on my_list'''
#     if my_list is None:
    #my_list=[0.1,0.2,0.5,1,2,5,10,20];

    L=len(my_list);
    x_label=np.zeros(L);
    mean_data=np.zeros(L);
    std_data=np.zeros(L);
    x0=0;
    for j,xj in enumerate(my_list):
        x_label[j]=(x0+xj)/2;
        index_0=(x>x0) & (x<=my_list[j]);

        # there are some satisfied elements;
        if sum(index_0)!=0:   
            x0=my_list[j];
            mean_data[j]=np.mean(y[index_0]);
            std_data[j]=np.std(y[index_0]);


    result=np.zeros((L,3));
    result[:,0]=x_label;
    result[:,1]=mean_data;
    result[:,2]=std_data;

    if plot_fig[0]:
        import my_fig_config_0
        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        ax.errorbar(result[:,0],result[:,1],result[:,2], fmt='.k',ecolor='gray')
        ax.set_xlabel("1T activity");
        if plot_fig[1]==0:
            ax.set_ylabel("2T expansion");
        elif plot_fig[1]==1:
            ax.set_ylabel("2T-R1 expansion");
        elif plot_fig[1]==2:
            ax.set_ylabel("2T-R2 expansion");

        ax.set_yscale('log')
        #ax.set_ylim((-0.1,10))
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')

        fig.tight_layout();

        if len(plot_fig)==1:
            fig.savefig(fig_dir+"1Tact_2Texp.eps")
        elif len(plot_fig)==2 & plot_fig[1]==1:
            fig.savefig(fig_dir+"1Tact_R1exp.eps")


    return result



def my_log_data(my_data):
    max_N=int(np.log(max(my_data))/np.log(2));
    log_data=np.zeros(max_N+1);
    for j in range(0,max_N+1):
        log_data[j]=np.sum((2**j<=my_data) & (my_data<2**(j+1)));
    
    return log_data;


def hamming(bc1,bc2): return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])

def match_or_not_0(bc1_long,bc2_long,N_HAMMING):
    #the method does assumes that the bc_list has been orderred, which may not be true in general
    # it is the updated version of match_or_not

    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    L1=len(bc1_list);
    L2=len(bc2_list);
    count_matrix=np.zeros((L1,L2));
    #if len(bc1_list)==len(bc2_list):
    for j in range(L1):
        for k in range(L2):
            if hamming(bc1_list[j],bc2_list[k])<=N_HAMMING:
                #print(count);
                count_matrix[j,k]=1;
    
    #the number of unique matches, in case there are redundant BCs in the list.  This redundancy is also removed in previous steps. So it is not necessary any more
    count=min(sum(count_matrix.sum(0)>0),sum(count_matrix.sum(1)>0))

    if count==L1 and count==L2:
        return True;
    else:
        return False;

def search_id(bc,bc_array):
    for j,bc_2 in enumerate(bc_array):
        if bc_2==bc: return j;
    return -1;

def match_or_not_1(bc1_long,bc2_long,N_HAMMING,dropout):
    #the method does assumes that the bc_list has been orderred, which may not be true in general
    # it is the updated version of match_or_not

    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    L1=len(bc1_list);
    L2=len(bc2_list);
    count_matrix=np.zeros((L1,L2));
    #if len(bc1_list)==len(bc2_list):
    for j in range(L1):
        for k in range(L2):
            if hamming(bc1_list[j],bc2_list[k])<=N_HAMMING:
                #print(count);
                count_matrix[j,k]=1;
    count=min(sum(count_matrix.sum(0)>0),sum(count_matrix.sum(1)>0))
    
    #considering the dropout, one cell may have less BC, but these BC should all match 

    xx=count/min(L1,L2); 
    yy=count/max(L1,L2); 
#    if xx>=1-0.00001 and yy>=1-dropout-0.00001:
    if yy>=1-dropout-0.00001:
#        print(count,L1,L2,xx,yy)
        # if yy>1: 
        #     print(bc1_long);
        #     print(bc2_long);
        return True;
    else:
        return False;

def match_N(bc1_long,bc2_long,N_HAMMING):
    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    L1=len(bc1_list);
    L2=len(bc2_list);
    count_matrix=np.zeros((L1,L2));
    #if len(bc1_list)==len(bc2_list):
    for j in range(L1):
        for k in range(L2):
            if hamming(bc1_list[j],bc2_list[k])<=N_HAMMING:
                #print(count);
                count_matrix[j,k]=1;
    count=min(sum(count_matrix.sum(0)>0),sum(count_matrix.sum(1)>0))
    
    #considering the dropout, one cell may have less BC, but these BC should all match 

    return count;

def match_or_not_4(bc1_long,bc2_long,N_HAMMING):
    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    L1=len(bc1_list);
    L2=len(bc2_list);
    count_matrix=np.zeros((L1,L2));
    #if len(bc1_list)==len(bc2_list):
    for j in range(L1):
        for k in range(L2):
            if hamming(bc1_list[j],bc2_list[k])<=N_HAMMING:
                #print(count);
                count_matrix[j,k]=1;
    count=min(sum(count_matrix.sum(0)>0),sum(count_matrix.sum(1)>0))
    
    #considering the dropout, one cell may have less BC, but these BC should all match 

    if count/max>=1:
        return True;
    else:
        return False;

def match_or_not_3(bc1_long,cell_N1,bc2_long,cell_N2,N_HAMMING,critical_cell_N,dropout):
    # add additional information of clone size 
    #considering the dropout, one cell may have less BC
    epsilon=0.0001;
    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    L1=len(bc1_list);
    L2=len(bc2_list);
    count_matrix=np.zeros((L1,L2));
    #if len(bc1_list)==len(bc2_list):
    for j in range(L1):
        for k in range(L2):
            if hamming(bc1_list[j],bc2_list[k])<=N_HAMMING:
                #print(count);
                count_matrix[j,k]=1;
    count=min(sum(count_matrix.sum(0)>0),sum(count_matrix.sum(1)>0))
 #   print(count)
#    count=int(count);
    

    # there is no dropout in both clones
    if cell_N2>=critical_cell_N  and cell_N1>=critical_cell_N  and count==L2 and count==L1:
        print("type 1:{},{},{}".format(count,L1,L2));
        return True;

    # there is no dropout in clone 1 but might be in clone 2 according to the clone size  
    if cell_N1>=critical_cell_N  and cell_N2<critical_cell_N  and count==L1 and count-L2>=0 and count-L2<=dropout: # allow at most 2 dropouts
        print("type 2:{},{},{}".format(count,L1,L2));
        return True;
    
    # there is no dropout in clone 2 but might be in clone 1 according to the clone size  
    if cell_N2>=critical_cell_N  and cell_N1<critical_cell_N  and count==L2 and count-L1>=0 and count-L1<=dropout: # allow at most 2 dropouts
        print("type 3:{},{},{}".format(count,L1,L2));
        return True;


    # there is dropout in both clones
    if cell_N2<critical_cell_N  and cell_N1<critical_cell_N:  
        if count==min(L1,L2) and max(L1,L2)-count<=dropout:
            print("type 4:{},{},{}".format(count,L1,L2));
            return True;
        else:
            return False;




def match_or_not_2(bc1_long,bc2_long,N_HAMMING,dropout):
    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    count_1=0;
    #if len(bc1_list)==len(bc2_list):
    for j in range(len(bc1_list)):
        for k in range(len(bc2_list)):
            if hamming(bc1_list[j],bc2_list[k])<=N_HAMMING:
                #print(count);
                count_1+=1;

    count_2=0;
    if len(bc1_list)==len(bc2_list):
        for j in range(len(bc1_list)):
            if hamming(bc1_list[j],bc2_list[j])<=N_HAMMING:
                #print(count);
                count_2+=1;
        
    xx=count_1/min(len(bc1_list),len(bc2_list)); #considering the dropout, one cell may have less BC, but these BC should all match 
    yy=count_1/max(len(bc1_list),len(bc2_list)); 
    if xx>=1-0.00001 and yy>=1-dropout-0.00001 and (count_2!=len(bc1_list) or count_2!=len(bc2_list)):
        print(count_1,count_2,len(bc1_list),len(bc2_list),xx,yy)
        print(bc1_long);
        print(bc2_long);
        return True;
    else:
        return False;

def match_or_not(bc1_long,bc2_long,N_HAMMING):
    bc2_list=bc2_long.split("_");
    bc1_list=bc1_long.split("_");
    count=0;

    #It is possible that the bc_list may contain several identical bc units.  
    # Here, we have a stringent requirement: the two clones has exactly the same array of bc units 
    # [each unit is similar enough within a certain Hamming distance] 

    #the method assumes that the bc_list has been orderred, which may not be true in general
    if len(bc1_list)==len(bc2_list):
        for j in range(len(bc1_list)):
            if hamming(bc1_list[j],bc2_list[j])<=N_HAMMING:
                #print(count);
                count+=1;
        
    if count==len(bc1_list):
        return True;
    else:
        return False;

def plot_figs_expansion_activity(cutoff_small,cutoff_large,pseudo_2,plot_fig,T1_HSC_info,T1_Kit_info,T2_HSC_info,T2_Kit_info,fig_dir):

    '''Here, you will need to provide the parameter: cutoff_large, T1_HSC_info, T1_Kit_info, T2_HSC_info,T2_Kit_info '''

#    pseudo_2=1;

    import my_fig_config_0

    used_clone_index=(T1_HSC_info>=cutoff_small)# & (T2_HSC_info>=cutoff_small);
    selected_T2_data=T2_HSC_info+T2_Kit_info;
    selected_T1_data=T1_HSC_info+T1_Kit_info;
    large_clone_index_T1=(selected_T1_data>cutoff_large) & used_clone_index
    large_clone_index_T2=(selected_T2_data>cutoff_large) & used_clone_index

    expansion_T1=(selected_T1_data+pseudo_2);
    potency_T1=(T1_Kit_info+pseudo_2)/(T1_HSC_info+pseudo_2);
    expansion_T2=(selected_T2_data+pseudo_2)/(T1_HSC_info+pseudo_2);
    potency_T2=(T2_Kit_info+pseudo_2)/(T2_HSC_info+pseudo_2);

    corr_T1actT2exp=pearsonr(potency_T1[used_clone_index],expansion_T2[used_clone_index])[0]
    corr_T1expT2exp=pearsonr(expansion_T1[used_clone_index],expansion_T2[used_clone_index])[0]
    corr_T1expT2act=pearsonr(expansion_T1[used_clone_index],potency_T2[used_clone_index])[0]
    corr_T1actT2act=pearsonr(potency_T1[used_clone_index],potency_T2[used_clone_index])[0]

    # ylabels
    if plot_fig[1]==0:
        my_ylabel_1="2T expansion";
        my_ylabel_2="2T activity";
    elif plot_fig[1]==1:
        my_ylabel_1="2T-R1 expansion";
        my_ylabel_2="2T-R1 activity";
    elif plot_fig[1]==2:
        my_ylabel_1="2T-R2 expansion";
        my_ylabel_2="2T-R2 activity";


    if plot_fig[0]==1:
        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        used_clone_index= (T1_HSC_info>=cutoff_small) 
        ax.plot(potency_T1[used_clone_index],expansion_T2[used_clone_index],".r")
        # ax.plot(potency_T1[used_clone_index & (T1_HSC_info==0)],expansion_T2[used_clone_index & (T1_HSC_info==0) ],".",color='gray')
        # #ax.plot(potency_T1[T2_HSC_info==0],expansion_T2[T1_HSC_info==0],".r")
        # ax.plot(potency_T1[large_clone_index_T1],expansion_T2[large_clone_index_T1],".b")
        # ax.plot(potency_T1[large_clone_index_T2],expansion_T2[large_clone_index_T2],".g")
        ax.set_xlabel("1T activity");
        if plot_fig[1]==0:
            ax.set_ylabel("2T expansion");
        if plot_fig[1]==1:
            ax.set_ylabel("2T-R1 expansion");
        if plot_fig[1]==2:
            ax.set_ylabel("2T-R2 expansion");
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        ax.text(0.65,0.85, 'C={:.2f}'.format(corr_T1actT2exp), fontsize=16,transform=ax.transAxes)
        fig.savefig(fig_dir+"potency_1T2T.eps")


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        used_clone_index= (T1_HSC_info>=cutoff_small) 
        ax.plot(expansion_T1[used_clone_index],expansion_T2[used_clone_index],".r")
        # ax.plot(expansion_T1[used_clone_index & (T1_HSC_info==0)],expansion_T2[used_clone_index & (T1_HSC_info==0)],".",color='gray')
        # ax.plot(expansion_T1[large_clone_index_T1],expansion_T2[large_clone_index_T1],".b")
        # ax.plot(expansion_T1[large_clone_index_T2],expansion_T2[large_clone_index_T2],".g")
        ax.set_xlabel("1T expansion");
        if plot_fig[1]==0:
            ax.set_ylabel("2T expansion");
        if plot_fig[1]==1:
            ax.set_ylabel("2T-R1 expansion");
        if plot_fig[1]==2:
            ax.set_ylabel("2T-R2 expansion");
        ax.set_yscale('log')
        ax.set_xscale('log')
        fig.tight_layout();
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        ax.text(0.65,0.85, 'C={:.2f}'.format(corr_T1expT2exp), fontsize=16,transform=ax.transAxes)
        fig.savefig(fig_dir+"expansion_1T2T.eps")


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        used_clone_index= (T2_HSC_info>=cutoff_small) 
        ax.plot(expansion_T1[used_clone_index],potency_T2[used_clone_index],".r")
        # ax.plot(expansion_T1[used_clone_index & (T2_HSC_info==0)],potency_T2[used_clone_index & (T2_HSC_info==0)],".",color='gray')
        # #ax.plot(expansion_T1[T1_HSC_info==0],potency_T2[T1_HSC_info==0],".",color='gray')
        # ax.plot(expansion_T1[large_clone_index_T1],potency_T2[large_clone_index_T1],".b")
        # ax.plot(expansion_T1[large_clone_index_T2],potency_T2[large_clone_index_T2],".g")
        ax.set_xlabel("1T expansion");
        if plot_fig[1]==0:
            ax.set_ylabel("2T activity");
        if plot_fig[1]==1:
            ax.set_ylabel("2T-R1 activity");
        if plot_fig[1]==2:
            ax.set_ylabel("2T-R2 activity");
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        ax.text(0.1,0.85, 'C={:.2f}'.format(corr_T1expT2act), fontsize=16,transform=ax.transAxes)
        fig.savefig(fig_dir+"expansion_activity_1T2T.eps")


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        used_clone_index=(T1_HSC_info>=cutoff_small) & (T2_HSC_info>=cutoff_small) 
        ax.plot(potency_T1[used_clone_index],potency_T2[used_clone_index],".r")
        # ax.plot(potency_T1[used_clone_index & (T1_HSC_info==0)],potency_T2[used_clone_index & (T1_HSC_info==0)],".",color='gray');
        # ax.plot(potency_T1[used_clone_index & (T2_HSC_info==0)],potency_T2[used_clone_index & (T2_HSC_info==0)],".",color='gray');
        # ax.plot(potency_T1[large_clone_index_T1],potency_T2[large_clone_index_T1],".b")
        # ax.plot(potency_T1[large_clone_index_T2],potency_T2[large_clone_index_T2],".g")
        ax.set_xlabel("1T activity");
        if plot_fig[1]==0:
            ax.set_ylabel("2T activity");
        if plot_fig[1]==1:
            ax.set_ylabel("2T-R1 activity");
        if plot_fig[1]==2:
            ax.set_ylabel("2T-R2 activity");
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        ax.text(0.1,0.85, 'C={:.2f}'.format(corr_T1actT2act), fontsize=16,transform=ax.transAxes)
        fig.savefig(fig_dir+"activity_1T2T.eps")

    result=np.zeros(4);
    result[0]=corr_T1actT2exp;
    result[1]=corr_T1expT2exp;
    result[2]=corr_T1expT2act;
    result[3]=corr_T1actT2act;
    return result


def plot_figs_correlation_normalized(plot_fig,pseudo_act_exp,T1_HSC_info,T1_Kit_info,T2_HSC_info,T2_Kit_info,fig_dir):
    '''Scatter plots of data in T1 and T2'''
    pseudo=0.1; # this number is just for plotting, not important to change

    T2_clone_info=T2_HSC_info+T2_Kit_info;
    T1_clone_info=T1_HSC_info+T1_Kit_info;
    norm_T2_HSC_info=(T2_HSC_info+pseudo_act_exp)/(T1_HSC_info+pseudo_act_exp);
    norm_T2_Kit_info=(T2_Kit_info+pseudo_act_exp)/(T1_HSC_info+pseudo_act_exp);
    norm_T2_clone_info=(T2_clone_info+pseudo_act_exp)/(T1_HSC_info+pseudo_act_exp);

    corr_T1HSC_T2HSC=pearsonr(T1_HSC_info,norm_T2_HSC_info)[0]
    corr_T1HSC_T2Kit=pearsonr(T1_HSC_info,norm_T2_Kit_info)[0]
    corr_T1Kit_T2Kit=pearsonr(T1_Kit_info,norm_T2_Kit_info)[0]
    corr_T1Kit_T2HSC=pearsonr(T1_Kit_info,norm_T2_HSC_info)[0]
    corr_T1HSC_T1Kit=pearsonr(T1_HSC_info,T1_Kit_info)[0]
    corr_T2HSC_T2Kit=pearsonr(norm_T2_HSC_info,norm_T2_Kit_info)[0]
    corr_T1cloneT2clone=pearsonr(T1_clone_info,norm_T2_clone_info)[0]

    import my_fig_config_0

    if plot_fig[0]==1:
        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        ax.plot(T1_HSC_info+pseudo,norm_T2_HSC_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_HSC_info[T1_HSC_info==0]+pseudo,norm_T2_HSC_info[T1_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_HSC_info[norm_T2_HSC_info==0]+pseudo,norm_T2_HSC_info[norm_T2_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        #ax.set_title("HSC -> HSC");
        if plot_fig[1]==0:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T HSC expansion");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R1 HSC expansion");
        if plot_fig[1]==2:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R2 HSC expansion");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1HSC_T2HSC), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T1-HSC-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_HSC_info+pseudo,norm_T2_Kit_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_HSC_info+pseudo,norm_T2_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_HSC_info[T1_HSC_info==0]+pseudo,norm_T2_Kit_info[T1_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_HSC_info[norm_T2_Kit_info==0]+pseudo,norm_T2_Kit_info[norm_T2_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T Kit expansion");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R1 Kit expansion");
        if plot_fig[1]==2:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R2 Kit expansion");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1HSC_T2Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T1-Kit-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,norm_T2_Kit_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_Kit_info+pseudo,norm_T2_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_Kit_info[T1_Kit_info==0]+pseudo,norm_T2_Kit_info[T1_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_Kit_info[norm_T2_Kit_info==0]+pseudo,norm_T2_Kit_info[norm_T2_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T Kit expansion");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R1 Kit expansion");
        if plot_fig[1]==2:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R2 Kit expansion");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1Kit_T2Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_Kit-T1-Kit-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,norm_T2_HSC_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_Kit_info+pseudo,norm_T2_HSC_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_Kit_info[T1_Kit_info==0]+pseudo,norm_T2_HSC_info[T1_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_Kit_info[norm_T2_HSC_info==0]+pseudo,norm_T2_HSC_info[norm_T2_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T HSC expansion");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R1 HSC expansion");
        if plot_fig[1]==2:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R2 HSC expansion");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1Kit_T2HSC), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_Kit-T1-HSC-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,norm_T2_HSC_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_HSC_info+pseudo,T1_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_HSC_info[T1_HSC_info==0]+pseudo,T1_Kit_info[T1_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_HSC_info[T1_Kit_info==0]+pseudo,T1_Kit_info[T1_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("1T Kit clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("1T Kit clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("1T Kit clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1HSC_T1Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T1-Kit-T1.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,norm_T2_HSC_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(norm_T2_HSC_info+pseudo,norm_T2_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(norm_T2_HSC_info[norm_T2_HSC_info==0]+pseudo,norm_T2_Kit_info[norm_T2_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(norm_T2_HSC_info[norm_T2_Kit_info==0]+pseudo,norm_T2_Kit_info[norm_T2_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("2T HSC expansion");
            ax.set_ylabel("2T Kit expansion");
        elif plot_fig[1]==1:
            ax.set_xlabel("R1 HSC expansion");
            ax.set_ylabel("R1 Kit expansion");
        if plot_fig[1]==2:
            ax.set_xlabel("R2 HSC expansion");
            ax.set_ylabel("R2 Kit expansion");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T2HSC_T2Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T2-Kit-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_clone_info+pseudo,norm_T2_clone_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_clone_info+pseudo,norm_T2_clone_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_clone_info[T1_clone_info==0]+pseudo,norm_T2_clone_info[T1_clone_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_clone_info[norm_T2_clone_info==0]+pseudo,norm_T2_clone_info[norm_T2_clone_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T clone size");
            ax.set_ylabel("2T expansion");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T clone size");
            ax.set_ylabel("2T-R1 expansion");
        if plot_fig[1]==2:
            ax.set_xlabel("1T clone size");
            ax.set_ylabel("2T-R2 expansion");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1cloneT2clone), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_T1-T2.eps");

    result=np.zeros(7);
    result[0]=corr_T1HSC_T2HSC;
    result[1]=corr_T1HSC_T2Kit
    result[2]=corr_T1Kit_T2Kit
    result[3]=corr_T1Kit_T2HSC
    result[4]=corr_T1HSC_T1Kit
    result[5]=corr_T2HSC_T2Kit
    result[6]=corr_T1cloneT2clone

    return result



def plot_figs_correlation(plot_fig,T1_HSC_info,T1_Kit_info,T2_HSC_info,T2_Kit_info,fig_dir):
    '''Scatter plots of data in T1 and T2'''
    pseudo=0.1; # this number is just for plotting, not important to change

    T2_clone_info=T2_HSC_info+T2_Kit_info;
    T1_clone_info=T1_HSC_info+T1_Kit_info;

    corr_T1HSC_T2HSC=pearsonr(T1_HSC_info,T2_HSC_info)[0]
    corr_T1HSC_T2Kit=pearsonr(T1_HSC_info,T2_Kit_info)[0]
    corr_T1Kit_T2Kit=pearsonr(T1_Kit_info,T2_Kit_info)[0]
    corr_T1Kit_T2HSC=pearsonr(T1_Kit_info,T2_HSC_info)[0]
    corr_T1HSC_T1Kit=pearsonr(T1_HSC_info,T1_Kit_info)[0]
    corr_T2HSC_T2Kit=pearsonr(T2_HSC_info,T2_Kit_info)[0]
    corr_T1cloneT2clone=pearsonr(T1_clone_info,T2_clone_info)[0]

    import my_fig_config_0

    if plot_fig[0]==1:
        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        ax.plot(T1_HSC_info+pseudo,T2_HSC_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_HSC_info[T1_HSC_info==0]+pseudo,T2_HSC_info[T1_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_HSC_info[T2_HSC_info==0]+pseudo,T2_HSC_info[T2_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        #ax.set_title("HSC -> HSC");
        if plot_fig[1]==0:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T HSC clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R1 HSC clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R2 HSC clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R1 HSC clone size");
            ax.set_ylabel("R2 HSC clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1HSC_T2HSC), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T1-HSC-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_HSC_info+pseudo,T2_Kit_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_HSC_info+pseudo,T2_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_HSC_info[T1_HSC_info==0]+pseudo,T2_Kit_info[T1_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_HSC_info[T2_Kit_info==0]+pseudo,T2_Kit_info[T2_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T Kit clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R1 Kit clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("2T-R2 Kit clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R1 HSC clone size");
            ax.set_ylabel("R2 Kit clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1HSC_T2Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T1-Kit-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,T2_Kit_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_Kit_info+pseudo,T2_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_Kit_info[T1_Kit_info==0]+pseudo,T2_Kit_info[T1_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_Kit_info[T2_Kit_info==0]+pseudo,T2_Kit_info[T2_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T Kit clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R1 Kit clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R2 Kit clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R1 Kit clone size");
            ax.set_ylabel("R2 Kit clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1Kit_T2Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_Kit-T1-Kit-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,T2_HSC_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_Kit_info+pseudo,T2_HSC_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_Kit_info[T1_Kit_info==0]+pseudo,T2_HSC_info[T1_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_Kit_info[T2_HSC_info==0]+pseudo,T2_HSC_info[T2_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T HSC clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R1 HSC clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T Kit clone size");
            ax.set_ylabel("2T-R2 HSC clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R1 Kit clone size");
            ax.set_ylabel("R2 HSC clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1Kit_T2HSC), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_Kit-T1-HSC-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,T2_HSC_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_HSC_info+pseudo,T1_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_HSC_info[T1_HSC_info==0]+pseudo,T1_Kit_info[T1_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_HSC_info[T1_Kit_info==0]+pseudo,T1_Kit_info[T1_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("1T Kit clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("1T Kit clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T HSC clone size");
            ax.set_ylabel("1T Kit clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R1 HSC clone size");
            ax.set_ylabel("R1 Kit clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1HSC_T1Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T1-Kit-T1.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_Kit_info+pseudo,T2_HSC_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T2_HSC_info+pseudo,T2_Kit_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T2_HSC_info[T2_HSC_info==0]+pseudo,T2_Kit_info[T2_HSC_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T2_HSC_info[T2_Kit_info==0]+pseudo,T2_Kit_info[T2_Kit_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("2T HSC clone size");
            ax.set_ylabel("2T Kit clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("R1 HSC clone size");
            ax.set_ylabel("R1 Kit clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("R2 HSC clone size");
            ax.set_ylabel("R2 Kit clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R2 HSC clone size");
            ax.set_ylabel("R2 Kit clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T2HSC_T2Kit), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_HSC-T2-Kit-T2.eps");


        fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
        #ax.plot(T1_clone_info+pseudo,T2_clone_info+pseudo,'.',color="r"); # detected only in T1
        ax.plot(T1_clone_info+pseudo,T2_clone_info+pseudo,'.',color="r"); # detected in both T1 and T2
        ax.plot(T1_clone_info[T1_clone_info==0]+pseudo,T2_clone_info[T1_clone_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        ax.plot(T1_clone_info[T2_clone_info==0]+pseudo,T2_clone_info[T2_clone_info==0]+pseudo,'.',color='gray'); # detected only in T1 or T2
        if plot_fig[1]==0:
            ax.set_xlabel("1T clone size");
            ax.set_ylabel("2T clone size");
        elif plot_fig[1]==1:
            ax.set_xlabel("1T clone size");
            ax.set_ylabel("2T-R1 clone size");
        if plot_fig[1]==2:
            ax.set_xlabel("1T clone size");
            ax.set_ylabel("2T-R2 clone size");
        if plot_fig[1]==3:
            ax.set_xlabel("R1 clone size");
            ax.set_ylabel("R2 clone size");
        ax.text(0.07,0.85, 'C={:.2f}'.format(corr_T1cloneT2clone), fontsize=16,transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if not_plot_upper_right_line:
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
        fig.tight_layout();
        fig.savefig(fig_dir+"clone_size_map_T1-T2.eps");

    result=np.zeros(7);
    result[0]=corr_T1HSC_T2HSC;
    result[1]=corr_T1HSC_T2Kit
    result[2]=corr_T1Kit_T2Kit
    result[3]=corr_T1Kit_T2HSC
    result[4]=corr_T1HSC_T1Kit
    result[5]=corr_T2HSC_T2Kit
    result[6]=corr_T1cloneT2clone

    return result



def correlation_scrambled_data(input_1,input_2):
    '''scramble the data to test the statistical significance'''
    tot_N=1000;
    corr=np.zeros(tot_N);
    for j in range(tot_N):
        np.random.shuffle(input_1);
        np.random.shuffle(input_2);
        corr[j]=pearsonr(input_1,input_2)[0];
    return corr




def inferred_clone_size(data,Pn):
    '''This is different from the negative binomial approach now'''
    data=data+1 # to convert to the negative binomial problem
    output=np.zeros(len(data),dtype="int")
    for j in range(len(data)):
        output[j]=data[j]+np.random.negative_binomial(data[j],Pn,1); #Pn the negative binomial success rate 
    output=output.astype(int);
    output=output-1; # convert back to our problem
    return output
    ## the negative binomial approach
    # pseudocount=1;
    # data=data+pseudocount;
    # output=np.zeros(len(data),dtype="int")
    # for j in range(len(data)):
    #     output[j]=data[j]+np.random.negative_binomial(data[j],Pn,1); #Pn the negative binomial success rate 
    # output=output.astype(int);
    # return output


def plot_histogram_with_a_line(simu_data,bin_N,exp_point,xlim,fig_dir,xlabel,file_name):
    
    Pvalue=sum(simu_data>exp_point)/len(simu_data);
    Pvalue=np.min([Pvalue,1-Pvalue]);
    print("Mean for "+xlabel+": ",np.mean(simu_data))
    print("Standard deviation for "+xlabel+": ",np.std(simu_data))
    print("Experimental value for "+xlabel+": ",exp_point)
    print("P value for "+xlabel+": ",Pvalue);


    import my_fig_config_0
    fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
    ax.hist(simu_data,bin_N,color="r");
    ax.set_xlabel(xlabel);
    ax.set_ylabel("Histogram");
    if len(xlim)==2:
        ax.set_xlim(xlim)

    xl = np.array(ax.get_xlim())
    yl = np.array(ax.get_ylim())
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.text(0.07,0.85, 'Pv={:.4f}'.format(Pvalue), fontsize=16,transform=ax.transAxes)
    ax.text(0.07,0.55, 'Exp={:.4f}'.format(exp_point), fontsize=16,transform=ax.transAxes)
    # Plot the threshold    
    ax.plot(exp_point* np.ones(2), yl, "--",c='b', linewidth=2)


    if not_plot_upper_right_line:
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

    fig.tight_layout();
    fig.savefig(fig_dir+file_name+".eps")



    return [np.mean(simu_data), np.std(simu_data),exp_point,Pvalue]



def plot_figs_1Tactivity_2Texpansion(cutoff_small,pseudo_2,threshold,plot_fig,T1_HSC_info,T1_Kit_info,T2_HSC_info,T2_Kit_info,fig_dir):
    used_clone_index=(T1_HSC_info>=cutoff_small)# & (T2_HSC_info>=cutoff_small);
    potency_T1=(T1_Kit_info+pseudo_2)/(T1_HSC_info+pseudo_2);
    expansion_T2=(T2_HSC_info+T2_Kit_info+pseudo_2)/(T1_HSC_info+pseudo_2);

    x1=potency_T1[used_clone_index];
    y1=expansion_T2[used_clone_index]

    data_1_y=y1[x1<=threshold]
    data_1_x=x1[x1<=threshold]
    data_2_y=y1[x1>threshold]
    data_2_x=x1[x1>threshold]


    import my_fig_config_0
    fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
    ax.plot(x1[x1>threshold],np.log(y1[x1>threshold]),'.',markersize=8,color='gray')
    ax.plot(x1[x1<=threshold],np.log(y1[x1<=threshold]),'.',markersize=8,color='gray')
    ax.set_xlabel("1T activity");
    if plot_fig[1]==0:
        ax.set_ylabel("2T expansion");
    elif plot_fig[1]==1:
        ax.set_ylabel("2T-R1 expansion");
    elif plot_fig[1]==2:
        ax.set_ylabel("2T-R2 expansion");
    #ax.set_yscale('log')
    #ax.set_ylim((-0.1,20))
    ax.set_xscale('log')
#   ax.set_yscale('log')
    xl = np.array(ax.get_xlim())
    yl = np.array(ax.get_ylim())
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    plt.yticks([np.log(0.01),np.log(0.1),np.log(1),np.log(10),np.log(100)], [r"$10^{-2}$", r"$10^{-1}$", r"$10^0$",r"$10^{1}$",r"$10^{2}$"])

    # Plot the threshold    
    ax.plot(threshold* np.ones(2), yl, "--",c='cyan', linewidth=2)
    #ax.plot(right_thresh* np.ones(2), yl,"--", c='black', linewidth=1)
    ax.errorbar(np.mean(data_1_x),np.log(np.mean(data_1_y)),np.std(np.log(data_1_y)), fmt='^r',markersize=7,ecolor='#FF499D')
    ax.errorbar(np.mean(data_2_x),np.log(np.mean(data_2_y)),np.std(np.log(data_2_y)), fmt='^r',markersize=7,ecolor='#FF499D') # #FFA69D

    if not_plot_upper_right_line:
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

    fig.tight_layout();
    fig.savefig(fig_dir+"potency_1T2T_mean_std_paper.eps")

    mean_ratio=np.mean(data_1_y)/np.mean(data_2_y);
    std_ratio=np.std(data_1_y)/np.std(data_2_y);

    print("data points # for second class:",len(data_2_x))
    print("mean_1/mean_2",mean_ratio)
    print("Std1/std2",std_ratio)


    #####
    my_list=[0.1,0.2,0.5,1.01,2,5,10,20];
    result=plot_average_expansion(x1,y1,my_list,plot_fig,fig_dir);

    return [mean_ratio,std_ratio]



def plot_figs_1Tactivity_2Texpansion_1(cutoff_small,pseudo_2,threshold,plot_fig,x1,y1,fig_dir):
    # used_clone_index=(T1_HSC_info>=cutoff_small)# & (T2_HSC_info>=cutoff_small);
    # potency_T1=(T1_Kit_info+pseudo_2)/(T1_HSC_info+pseudo_2);
    # expansion_T2=(T2_HSC_info+T2_Kit_info+pseudo_2)/(T1_HSC_info+pseudo_2);

    # x1=potency_T1[used_clone_index];
    # y1=expansion_T2[used_clone_index]
    correct_fig=1;

    data_1_y=y1[x1<=threshold]
    data_1_x=x1[x1<=threshold]
    data_2_y=y1[x1>threshold]
    data_2_x=x1[x1>threshold]


    import my_fig_config_0
    fig = plt.figure(figsize=(4,3.5)); ax = fig.add_subplot(1, 1, 1);
    if correct_fig:
        ax.plot(x1[x1>threshold],np.log(y1[x1>threshold]),'.',markersize=8,color='gray')
        ax.plot(x1[x1<=threshold],np.log(y1[x1<=threshold]),'.',markersize=8,color='gray')
    else:
        ax.plot(x1[x1>threshold],y1[x1>threshold],'.',markersize=8,color='gray')
        ax.plot(x1[x1<=threshold],y1[x1<=threshold],'.',markersize=8,color='gray')
    ax.set_xlabel("1T activity");
    if plot_fig[1]==0:
        ax.set_ylabel("2T expansion");
    elif plot_fig[1]==1:
        ax.set_ylabel("2T-R1 expansion");
    elif plot_fig[1]==2:
        ax.set_ylabel("2T-R2 expansion");
#    ax.set_yscale('log')
    #ax.set_ylim((-0.1,20))
    ax.set_xscale('log')
    if correct_fig != 1:
       ax.set_yscale('log')
    xl = np.array(ax.get_xlim())
    yl = np.array(ax.get_ylim())
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    if correct_fig:
        plt.yticks([np.log(0.01),np.log(0.1),np.log(1),np.log(10),np.log(100)], [r"$10^{-2}$", r"$10^{-1}$", r"$10^0$",r"$10^{1}$",r"$10^{2}$"])

    # Plot the threshold    
    ax.plot(threshold* np.ones(2), yl, "--",c='cyan', linewidth=2)
    #ax.plot(right_thresh* np.ones(2), yl,"--", c='black', linewidth=1)
    ax.errorbar(np.mean(data_1_x),np.log(np.mean(data_1_y)),np.std(np.log(data_1_y)), fmt='^r',markersize=7,ecolor='#FF499D')
    ax.errorbar(np.mean(data_2_x),np.log(np.mean(data_2_y)),np.std(np.log(data_2_y)), fmt='^r',markersize=7,ecolor='#FF499D') # #FFA69D

    if not_plot_upper_right_line:
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

    fig.tight_layout();
    fig.savefig(fig_dir+"potency_1T2T_mean_std_paper.eps")

    mean_ratio=np.mean(data_1_y)/np.mean(data_2_y);
    std_ratio=np.std(data_1_y)/np.std(data_2_y);

    print("data points # for second class:",len(data_2_x))
    print("mean_1/mean_2",mean_ratio)
    print("Std1/std2",std_ratio)


    #####
    my_list=[0.1,0.2,0.5,1.01,2,5,10,20];
    result=plot_average_expansion(x1,y1,my_list,plot_fig,fig_dir);

    return [mean_ratio,std_ratio]



