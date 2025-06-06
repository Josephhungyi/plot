#!/usr/bin/python3

import sys
import pdfplumber 
import pandas as pd
import numpy as np
import csv
import os
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import shutil
from datetime import datetime
import math

"""
class CMT(base,ChrMode,CheckRef): 
    category="exon,chromosome"
    plot="ExonName,CheckRef,ChrMode,TimeExon"
    arrow="yes,no"
    PeakSetChange="yes,no"
    BarSetChange="yes,no"
    DataOrderChange="yes,no"
    def classification(self, element,sample_name):
    def merge(self,upload_folder):
    def peak_set(self,sample_name,ax,color_list):
    def bar_set(self,sample_name,ax,color_list):
    def data_order_change(self,df): return df
"""
#分類程式
class base:
    arrow="yes"
    PeakSetChange="no"
    BarSetChange="no"
    DataOrderChange="no"
    
class ChrMode:
    def classification(self, element,sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if i==0:
                grp_check1=element[2][:2]
                grp_check2=element[2][:2]
            if remainder==2:
                grp_check2=element[2][:2]
                if grp_check1!=grp_check2 and 'Xq' not in grp_check2 and 'Yq' not in grp_check2:
                    grp_num+=1
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write(str(grp_num)+','+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write(str(grp_num)+','+element[j]+',')
                grp_check1=grp_check2
            elif remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')            
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')

class CheckRef:
    def classification(self, element,sample_name):
        for j in range(len(element)): 
            remainder=j%12
            # print(remainder,)
            if remainder==11:
                # print(element[j])
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('CHEK'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('check,'+element[j]+',')            
                elif element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('reference,'+element[j]+',')
                elif element[j].startswith('BRCA1'):
                    exon_num=element[j].split('-')[1]
                    if exon_num.isnumeric() and int(exon_num)>=5 and int(exon_num)<=24:
                        with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                            f.write('sample,BRCA1-'+str(int(exon_num)-1)+',')
                    else:
                         with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                            f.write('sample,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('sample,'+element[j]+',')                        
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')
        


#各分析,照字母排列

class ADPKD(base):
    category="exon"
    plot="ExonName"
    def classification(self, element,sample_name): 
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1: 
                if element[j].startswith('TSC2'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('PKD1'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('2,'+element[j]+',')
                elif element[j].startswith('PKD2'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('3,'+element[j]+',')
                elif element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('5,'+element[j]+',') 
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')
    def merge(self,upload_folder):
        txt_files=os.listdir(upload_folder+'_txt')
        txt_files.sort()
        for i in range(0,len(txt_files),2):
            sample_n=txt_files[i].split('-SET')[0]
            df1 = pd.read_csv(upload_folder+'_txt/'+sample_n+'-SET1.txt',header=0)
            
            df2 = pd.read_csv(upload_folder+'_txt/'+sample_n+'-SET2.txt',header=0)
            
            df=pd.concat([df1,df2])
            new_df=df.groupby(df.Category)
            TSC2_df=new_df.get_group(1).sort_values(by='hg18_loc.')
            PKD1_df=new_df.get_group(2).sort_values(by='hg18_loc.')
            PKD2_df=new_df.get_group(3).sort_values(by='hg18_loc.')
            ref_df=new_df.get_group(5).sort_values(by='hg18_loc.')
            merge_df=pd.concat([TSC2_df,PKD1_df,PKD2_df,ref_df])
            merge_txt=upload_folder+'_txt/'+sample_n+'.txt'
            merge_df.to_csv(merge_txt)
            barplot(sample_n,mode)
class ANUE(base,ChrMode):
    category="chromosome"
    plot="CheckRef"

class ATM(base,CheckRef):
    category="exon"
    plot="CheckRef"
    def merge(self,upload_folder):
        txt_files=os.listdir(upload_folder+'_txt')
        txt_files.sort()
        for i in range(0,len(txt_files),2):
            sample_n=txt_files[i].split('-SET')[0]
            # print(sample_n)
            df1 = pd.read_csv(upload_folder+'_txt/'+sample_n+'-SET1.txt',header=0)
           
            df2 = pd.read_csv(upload_folder+'_txt/'+sample_n+'-SET2.txt',header=0)
            
            df=pd.concat([df1,df2])
            new_df=df.groupby(df.Category)     
            
            new_sample=new_df.get_group('sample').sort_values(by='hg18_loc.')
            new_ref=new_df.get_group('reference').sort_values(by='hg18_loc.')
            merge_df=pd.concat([new_sample,new_ref])
            merge_txt=upload_folder+'_txt/'+sample_n+'.txt'
            merge_df.to_csv(merge_txt)
            # print(new_sample,new_ref)
            barplot(sample_n,mode)
            
    
class BRCA(base,CheckRef):
    category="exon"
    plot="CheckRef" 

class CATCH(base,ChrMode):
    category="chromosome"
    plot="ChrMode"
class CMT(base): 
    category="exon"
    plot="ExonName"
    def classification(self, element,sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('5,'+element[j]+',')            
                elif element[j].startswith('COX') or element[j].startswith('PPP2R3B'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('PMP'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('2,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                            f.write('3,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')

class MECP(base):
    category="exon"
    plot="ExonName"
    def classification(self, element,sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('5,'+element[j]+',')            
                elif element[j].startswith('NTNG'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('CDKL'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('2,'+element[j]+',')
                elif element[j].startswith('ARX'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('4,'+element[j]+',')                        
                elif element[j].startswith('MECP2'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('6,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                            f.write('3,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')


class MLH(base):
    category="exon"
    plot="ExonName"
    def classification(self, element,sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('4,'+element[j]+',')            
                elif element[j].startswith('EPCAM'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('MSH'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('2,'+element[j]+',')
                elif element[j].startswith('MLH'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('3,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('4,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')

class MSH(base):
    category="exon"
    plot="ExonName"
    def classification(self, element, sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('5,'+element[j]+',')            
                elif element[j].startswith('MUTYH'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('EPCAM'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('2,'+element[j]+',')
                elif element[j].startswith('MLH'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('4,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                            f.write('3,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')

class NF1(base,CheckRef):
    category="exon"
    plot="CheckRef"
    
    
class PMS(base): 
    category="exon"
    plot="ExonName"
    def classification(self, element,sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('5,'+element[j]+',')            
                elif element[j].startswith('PMS'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('3,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')
        
 

class SHOX(base):
    category="exon"
    plot="ExonName"
    def classification(self, element, sample_name):
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:
                if element[j].startswith('Reference'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('5,'+element[j]+',')            
                elif element[j].startswith('SHOX-AREA') or element[j].startswith('PPP2R3B'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('SHOX'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('2,'+element[j]+',')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('3,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')


class SMN(base):
    category="exon"
    plot="ExonName" #TimeExon,ExonName
    arrow="no"
    PeakSetChange="yes"
    BarSetChange="yes"
    DataOrderChange="yes"
    def classification(self, element,sample_name): 
        for j in range(len(element)): 
            remainder=j%12
            if remainder==11:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+'\n')
            elif remainder==1:                                
                if element[j].startswith('SMN1&2-Int') or element[j].startswith('SMN1&2-8'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('1,'+element[j]+',')
                elif element[j].startswith('SMN1-7') or element[j].startswith('SMN2-7'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('9,'+element[j][0:5]+'E'+element[j][5]+',')
                elif element[j].startswith('SMN2-8') or element[j].startswith('SMN1-8'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('8,'+element[j][0:5]+'E'+element[j][5]+',')        
                elif element[j].startswith('NAIP'):
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                        f.write('7,Reference,')
                else:
                    with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                            f.write('7,'+element[j]+',')
            else:
                with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                    f.write(element[j]+',')
    def peak_set(self,sample_name,ax,color_list):
        ax.set_title(sample_name,pad=10,loc='center',fontsize=12)
        custom_lines = [Line2D([0], [0], color=color_list[1], lw=4),
                Line2D([0], [0], color=color_list[9], lw=4),
                Line2D([0], [0], color=color_list[8], lw=4)]
        plt.subplots_adjust(top=0.95)
        plt.legend(custom_lines, ['SMN1+SMN2', 'Exon7', 'Exon8'], loc='upper right',fontsize=12)
        plt.savefig(results_folder+'/'+sample_name+'_2.jpg')
    def bar_set(self,sample_name,ax,color_list):
        ax.set_title(sample_name,pad=10,loc='center',fontsize=12)
        custom_lines = [Line2D([0], [0], color=color_list[1], lw=4),
                Line2D([0], [0], color=color_list[9], lw=4),
                Line2D([0], [0], color=color_list[8], lw=4)]
        plt.subplots_adjust(top=0.92)
        plt.legend(custom_lines, ['SMN1+SMN2', 'Exon7', 'Exon8'], loc='upper right',fontsize=12)
        plt.savefig(results_folder+'/'+sample_name+'_3.jpg')
    def data_order_change(self,df):
        reorder_date=df.iloc[0:4]
        reorder_date=reorder_date.sort_values(by='Dnt')
        reorder_date=reorder_date.reset_index(drop=True)
        df.iloc[0:4]=reorder_date
        return df
        


def pdf2table(pid,file,mode):
   
    chr_list=['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','Xp','Xq','Yp','Yq']
    
    selected_class=getattr(sys.modules[__name__], mode)
    
    if file.endswith('.pdf'):
        # print(i)
        sample_name=file.split('.')[0]
    #pdf轉txt
    pdf = pdfplumber.open(upload_folder+'/'+sample_name+'.pdf')
    page0 = pdf.pages[0]
    header=page0.extract_table()[2][1].replace('[','').replace(']','').split()
    header=header[1:]
    gender=page0.extract_table()[2][0].split('gender: ')[1].split(' RPQ:')[0]
    table = page0.extract_table()[3][1].replace('%','')
    table = page0.extract_table()[3][1].replace('%','').split('\n')
    table = table[:-1]   
    n=0
    for i in header:
        if i.startswith('nt'):
            header[n-1:n+1] = [''.join(header[n-1:n+1])]
        elif i.startswith('loc'):
            header[n-1:n+1] = ['_'.join(header[n-1:n+1])]
        elif i.startswith('detail'):
            header[n-1:n+1] = ['_'.join(header[n-1:n+1])]
        n+=1   
    if selected_class.category=="exon":
        header.insert(1,'Category')
    elif selected_class.category=="chromosome":
        header.insert(2,'Category')
    file_exsit=os.path.exists(upload_path+'/'+folder_name+'_txt/'+sample_name+'.txt')
    if file_exsit is True:
        os.remove(upload_folder+'_txt/'+sample_name+'.txt')
    for i in header:
        if i == header[-1]:
            with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                f.write(i+'\n')
        else:            
            with open(upload_folder+'_txt/'+sample_name+'.txt','a+') as f:
                f.write(i+',')

    grp_num=0
    #table整理
    for i in range(len(table)):
        element=table[i].split()
        if len(element)==12 and element[1].endswith('Int') and element[2].isdigit:              
            element[1:3] = ['_'.join(element[1:3])]
            element.insert(7,str(0))
        elif len(element)==11 and len(element[1]) > 14:
            element.insert(2,element[1][14:])
            element[1]=element[1][:14]
        elif len(element)==11 and len(element[1]) <= 14:
            element.insert(7,str(0))
            # print(element[1],element[2])
        elif len(element)==13:
            if element[2][0:2] not in chr_list:
                element[1:3] = ['_'.join(element[1:3])]
            else:
                element[11:13] = ['_'.join(element[11:13])]
        elif len(element)==14:
            element[1:3] = ['_'.join(element[1:3])]
            element[11:13] = ['_'.join(element[11:13])]
                        
        #依模式分類
        selected_class().classification(element,sample_name)
                   
    os.chmod(upload_folder+'_txt/'+sample_name+'.txt', 0o777)
    return sample_name,gender  
    
    
    
def barplot(sample_name,mode):
    selected_class=getattr(sys.modules[__name__], mode)
    #color_list=['#73edde','#9966ff','#6699ff','#f490e8','#ffd966','#7df383','#ff7c80','#afabab','#04AC00']
    color_list=['#73edde','#9966ff','#6699ff','#f490e8','#ffd966','#7df383','#ff7c80','#afabab','#03ac00','#17cbff','#04AC00']
    # 畫圖            
    df = pd.read_csv(upload_folder+'_txt/'+sample_name+'.txt',header=0)

    df['Height'] = pd.to_numeric(df['Height'])
    df['Dnt'] = pd.to_numeric(df['Dnt'])

    #bar圖
    plt.figure(figsize=(15,6))
    plt.title(sample_name,fontsize=25)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('gray')
    ax.spines['top'].set_color('none')
    ax.tick_params(length=0)
    #背景色
    # ax.set_facecolor('#AEAED6') 
    #輔助線
    plt.grid(axis='y',c='gray',zorder=0)
    #下緣寬度設定
    plt.gcf().subplots_adjust(bottom=0.25)
    #x軸文字轉角度
    # plt.setp( ax.xaxis.get_majorticklabels(), rotation=60, ha="right" )
    plt.xticks(rotation=90)
    plt.ylabel('Ratio',fontsize=14)
    plt.ylim([0,2])
    #間距刻度
    plt.yticks(np.arange(0,2.5,0.5))
    
    if selected_class.plot == "CheckRef":
        barlist=plt.bar(np.arange(len(df.index)),df.loc[:,'Ratio'],width=0.8,linewidth=1,tick_label=df.loc[:,'Gene-Exon'], zorder=3)
        for i in range(len(df.index)):
            if df.loc[i,'Category']=='sample':
                barlist[i].set_color('#73edde')
                # barlist[i].set_edgecolor('black')
            elif df.loc[i,'Category']=='check':
                barlist[i].set_color('#9966ff')
                # barlist[i].set_edgecolor('black')
            elif df.loc[i,'Category']=='reference':
                barlist[i].set_color('#6699ff')
                # barlist[i].set_edgecolor('black')
        plt.margins(0)
    elif selected_class.plot == "ChrMode":
        label_list=[]
        for i in df.loc[:,'Chr.band']:
            if i.startswith('0'):
                label_list.append('Chromosome '+i[1])
            elif i.startswith('X'):
                label_list.append('Chromosome X')
            elif i.startswith('Y'):
                label_list.append('Chromosome Y')
            else:
                label_list.append('Chromosome '+i[0:2])
        barlist=plt.bar(np.arange(len(df.index)),df.loc[:,'Ratio'],width=0.8,linewidth=1,tick_label=label_list, zorder=3)
        for i in range(len(df.index)):
            color_num=df.loc[i,'Category']
            barlist[i].set_color(color_list[int(color_num)])
            
        
    elif selected_class.plot == "ExonName":
        if selected_class.DataOrderChange == "yes":
            df=selected_class().data_order_change(df)
            
        barlist=plt.bar(np.arange(len(df.index)),df.loc[:,'Ratio'],width=0.8,linewidth=1,tick_label=df.loc[:,'Gene-Exon'], zorder=3)
        for i in range(len(df.index)):
            color_num=df.loc[i,'Category']
            barlist[i].set_color(color_list[int(color_num)])
            
    elif selected_class.plot == "TimeExon":
        df=df.sort_values(by='Dnt')
        df=df.reset_index(drop=True)
        barlist=plt.bar(np.arange(len(df.index)),df.loc[:,'Ratio'],width=0.8,linewidth=1,tick_label=df.loc[:,'Gene-Exon'], zorder=3)
        for i in range(len(df.index)):
            color_num=df.loc[i,'Category']
            #print(color_num,df.loc[i,'Gene-Exon'])
            barlist[i].set_color(color_list[int(color_num)])
    # plt.margins(0)
    plt.margins(0.005)
    if selected_class.BarSetChange == 'yes':
        selected_class().bar_set(sample_name,ax,color_list)
    else:
        plt.savefig(results_folder+'/'+sample_name+'_bar.jpg')
    plt.close()

    # plt.figure(figsize=(10,6))
def peakplot(sample_name,mode,gender):
    selected_class=getattr(sys.modules[__name__], mode)
    color_list=['#73edde','#9966ff','#6699ff','#f490e8','#ffd966','#7df383','#ff7c80','#afabab','#03ac00','#17cbff']
    df = pd.read_csv(upload_folder+'_txt/'+sample_name+'.txt',header=0)
    df['Height'] = pd.to_numeric(df['Height'])
    df['Dnt'] = pd.to_numeric(df['Dnt'])
   
    #peak圖
    max_value=df['Height'].max()
    if mode=="SMN":
        fig, ax = plt.subplots(figsize=(12,8), dpi= 80)
    else:
        fig, ax = plt.subplots(figsize=(15,8), dpi= 80)
    for i in range(len(df.index)):
        string_len=len(df.loc[i,'Gene-Exon'])
        ratio=math.log(270/int(string_len))
        #畫箭頭
        if selected_class.arrow=='no':
            pass
        elif 'Y' in df.loc[i,'Chr.band'] and gender=='Female':
            if df.loc[i,'Ratio'] >0:
                if len(df.loc[i,'Gene-Exon']) > 7:
                    plt.arrow(df.Dnt[i],df.Height[i]+max_value/ratio,0,-(max_value/20),width=1.5,color="r",head_length=max_value/40)
                else:
                    plt.arrow(df.Dnt[i],df.Height[i]+max_value/4.3,0,-(max_value/20),width=1.5,color="r",head_length=max_value/40)
        elif df.loc[i,'Gene-Exon'].startswith('MSH2-10Mbinv'):
            if df.loc[i,'Ratio'] >0:
                if len(df.loc[i,'Gene-Exon']) > 7:
                    plt.arrow(df.Dnt[i],df.Height[i]+max_value/ratio,0,-(max_value/20),width=1.5,color="r",head_length=max_value/40)
                else:
                    plt.arrow(df.Dnt[i],df.Height[i]+max_value/4.3,0,-(max_value/20),width=1.5,color="r",head_length=max_value/40)
        elif 'MUT' not in df.loc[i,'Gene-Exon']:
            if df.loc[i,'Ratio'] <= 0.65 or df.loc[i,'Ratio'] >= 1.25:
                if len(df.loc[i,'Gene-Exon']) > 7:
                    plt.arrow(df.Dnt[i],df.Height[i]+max_value/ratio,0,-(max_value/20),width=1.5,color="r",head_length=max_value/40)
                else:
                    plt.arrow(df.Dnt[i],df.Height[i]+max_value/4.3,0,-(max_value/20),width=1.5,color="r",head_length=max_value/40)

        #畫peak
        if selected_class.plot == "CheckRef":
            plt.text(df.Dnt[i]-1.5,df.Height[i]+max_value/40,df.loc[i,'Gene-Exon'],rotation=90)            
            if df.loc[i,'Category']=='sample':
                ax.vlines(x=df.Dnt[i], ymin=0, ymax=df.Height[i], color='#73edde', linewidth=2)
                ax.scatter(x=df.Dnt[i], y=df.Height[i], s=75, color='#73edde')
            elif df.loc[i,'Category']=='check':
                ax.vlines(x=df.Dnt[i], ymin=0, ymax=df.Height[i], color='#9966ff', linewidth=2)
                ax.scatter(x=df.Dnt[i], y=df.Height[i], s=75, color='#9966ff')
            elif df.loc[i,'Category']=='reference':
                ax.vlines(x=df.Dnt[i], ymin=0, ymax=df.Height[i], color='#6699ff', linewidth=2)
                ax.scatter(x=df.Dnt[i], y=df.Height[i], s=75, color='#6699ff')
        elif selected_class.plot == "ChrMode":
            df['Category'] = pd.to_numeric(df['Category'])
            plt.text(df.Dnt[i]-1.5,df.Height[i]+max_value/40,df.loc[i,'Chr.band'],rotation=90)            
            color_num=df.loc[i,'Category']
            ax.vlines(x=df.Dnt[i], ymin=0, ymax=df.Height[i], color=color_list[int(color_num)], linewidth=2)
            ax.scatter(x=df.Dnt[i], y=df.Height[i], s=75, color=color_list[int(color_num)])
        elif selected_class.plot == "ExonName" or selected_class.plot == "TimeExon":
            df['Category'] = pd.to_numeric(df['Category'])
            plt.text(df.Dnt[i]-1.5,df.Height[i]+max_value/40,df.loc[i,'Gene-Exon'],rotation=90)            
            color_num=df.loc[i,'Category']
            ax.vlines(x=df.Dnt[i], ymin=0, ymax=df.Height[i], color=color_list[int(color_num)], linewidth=2)
            ax.scatter(x=df.Dnt[i], y=df.Height[i], s=75, color=color_list[int(color_num)])
   
    ax.set_ylabel('RFU')
    ax.set_xticks(df.Dnt)
    plt.xticks(rotation=90)
    # plt.gca().set_ylim(bottom=0)
    # plt.gca().set_ylim(bottom=0)
    plt.ylim([0,max_value+max_value/2.5])   
    plt.yticks(np.arange(0,max_value+max_value/2.5,1000))
    # print(max_value)
    #版面設置
    if selected_class.PeakSetChange == 'yes':
        selected_class().peak_set(sample_name,ax,color_list)
    else:
        ax.set_title(sample_name, y=1.0, pad=-14, loc='right')
        plt.savefig(results_folder+'/'+sample_name+'_peak.jpg')
 
    plt.close()

if __name__ == '__main__':
    pid=sys.argv[1] 
    input_zip=sys.argv[2] 
    mode=sys.argv[3]
    merged=sys.argv[4]
    selected_class=getattr(sys.modules[__name__], mode)
    nowtime = datetime.now().strftime("%Y%m%d")
    upload_path="/the/path/"
    results_path="/the/path"
    folder_name=input_zip.split('.')[0]

    upload_folder=upload_path+'/'+folder_name
    outname=nowtime+'_'+mode+'_'+folder_name
    results_folder=results_path+'/'+outname
    check_folder=results_path+'/'+nowtime+'_'+folder_name
    n=1
    file_exsit=os.path.isdir(results_folder)
    if file_exsit is True:
        shutil.rmtree(results_folder, ignore_errors=True)
    # while os.path.exists(results_folder):
        # check_folder='_'.join([results_folder,str(n)])
        # n+=1
    # results_folder=check_folder
    os.mkdir(results_folder)
    if not os.path.exists(upload_folder+'_txt'):
        os.mkdir(upload_folder+'_txt')   
    files=os.listdir(upload_folder)
    
    #正常出圖
    #for i in files:
    #    [sample_name,gender]=pdf2table(pid,i,mode)
    #    peakplot(sample_name,mode,gender)
    #    barplot(sample_name,mode)
    #if merged=='Yes':
    #    selected_class().merge(upload_folder)
    
    #SMN不出bar圖
    for i in files:
        [sample_name,gender]=pdf2table(pid,i,mode)
        peakplot(sample_name,mode,gender)
        if merged=='No' and mode == 'SMN':
            pass
        elif merged=='No' and mode != 'SMN':       
            barplot(sample_name,mode)
    if merged=='Yes':
        selected_class().merge(upload_folder)
            
    shutil.make_archive(results_folder, 'zip', base_dir=outname,root_dir=results_path)
