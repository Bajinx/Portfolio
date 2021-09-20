#%%
import pandas
import pyBigWig
import numpy as np
import sys
import random
import pybedtools
import subprocess
import seaborn as sns
from matplotlib import pyplot as plt
from datetime import datetime



#%% Welche gene an welchen TF binden
# Read files and build TF database
promotor_bed_file_location = "/home/chris/Desktop/Ressources/Human_genome_files/Promotors_from_ensemble/Promotors_from_ensemble_most_5_no_pseudo.bed"  # sys.argv[2]
tf_db = pandas.read_csv("/home/chris/Desktop/Ressources/TF_Chipseqs/TF_DB.csv")
tf_db["Cellline"]=[str(i).strip() for i in tf_db["Cellline"]]
list_of_all_tfs=tf_db[tf_db["Cellline"]=="K562"]["TF"].tolist()

#%% Finding regions that are bound
dict_of_tf_with_bound_genes={}
for tf in list_of_all_tfs[:]:
    acc = tf_db[tf_db["TF"] == tf]["Peaks"].to_list()[0].strip() # Access identifier for that TF
    print(acc)
    peak_bed_file = "/home/chris/Desktop/Ressources/TF_Chipseqs/Peak_calling_files/" + acc + ".bed"
    peak_bed_file_df = pandas.read_csv("/home/chris/Desktop/Ressources/TF_Chipseqs/Peak_calling_files/" + acc + ".bed",delimiter="\t",header=None)
    bed_file_of_all_genes_df = pandas.read_csv(promotor_bed_file_location, "\t", header=None)
    bed_file_of_all_genes_df.columns = ["chr", "start", "end", "gene", "score", "strand"]

    #Overlap of Bed file with promoter regions an TF regions
    command = 'bedtools intersect -u -a ' + promotor_bed_file_location + ' -f 1.0 -F 1.0 -e -b ' + peak_bed_file
    a = subprocess.check_output(command.split()).decode("UTF-8")
    a = a.split("\n")
    b=pandas.DataFrame([x.split('\t') for x in a])
    b.dropna(inplace=True)
    dict_of_tf_with_bound_genes[tf]=b.iloc[:,3] #Resulting dict: Key TF; value all bound promoter names


#%%
import pickle

pickle.dump(initial_genesets_going_into_analysis, open("/home/chris/Desktop/Ressources/Analysen/From_chip_to_genelist/Analyse_all_ensemble_promotoren/initial_genesets_going_into_analysis_5000.p", "wb"))

#%% How many promoters are there for each TF (list and as plot)
list_of_gene_set_sizes=[]
counter=0
for key,value in dict_of_tf_with_bound_genes.items():
    if counter<3:
        list_of_gene_set_sizes.append(value.shape[0])

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
sns.swarmplot(list_of_gene_set_sizes,color="black",ax=ax1,c = [30 for i in range(len(list_of_gene_set_sizes))])
ax1.axvline(x=200,color="grey")
ax1.set_xlabel("Number of genes in the genesets")
#ax1.set_ylim([-1.7,1.7])
fig.show()


#%% only those TFs that have more than 200 promoter hits
dict_of_tf_with_bound_genes_200={}
for key,value in dict_of_tf_with_bound_genes.items():
    if 200<len(value.tolist()):
        print(key,value.tolist())
        dict_of_tf_with_bound_genes_200[key]=value.tolist()

#%% All TF reformating
dict_of_tf_with_bound_genes_all={}
for key,value in dict_of_tf_with_bound_genes.items():
        print(key,value.tolist())
        dict_of_tf_with_bound_genes_all[key]=value.tolist()

#%%
import subprocess

for initial_gene in list_of_initial_genes:

        print(initial_gene)
        # Motif and Database locations for all tools
        location="/home/chris/Desktop/Ressources/Analysen/From_chip_to_genelist/Analyse_all_ensemble_promotoren/19_02_fewer_tool_variants_2_ranking_variants/"+initial_gene
        tf_db = pandas.read_csv("/home/chris/Desktop/Ressources/TF_Chipseqs/TF_DB.csv")

        set_of_transfac_Motifs="/home/chris/Desktop/Ressources/Motifs/Hocomoco_Human/Core/HOCOMOCOv11_core_HUMAN_mono_transfac_format_filtered" #set of motifs in transfac format (JASPAR2020_HUMAN_transfac_P0.txt)
        set_of_meme_Motifs="/home/chris/Desktop/Ressources/Motifs/Hocomoco_Human/Core/HOCOMOCOv11_core_HUMAN_mono_meme_format_filtered.meme"
        #set_of_homer_Motifs_001="/home/chris/Desktop/Ressources/Motifs/Hocomoco_Human/Core/HOMER/HOCOMOCOv11_core_HUMAN_mono_homer_format_0.001.motif"
        #set_of_homer_Motifs_0005="/home/chris/Desktop/Ressources/Motifs/Hocomoco_Human/Core/HOMER/HOCOMOCOv11_core_HUMAN_mono_homer_format_0.0005.motif"
        set_of_homer_Motifs_0001="/home/chris/Desktop/Ressources/Motifs/Hocomoco_Human/Core/HOMER/HOCOMOCOv11_core_HUMAN_mono_homer_format_0.0001_filtered.motif"
        #path_to_EPIREGIO_Repo="/home/chris/Desktop/EPIREGIO/PASTAA_26_06_2020/ApplicationScenarioExamples/" #path to cloned gitHub repo
        set_of_jaspar_Motifs="/home/chris/Desktop/Ressources/Motifs/Hocomoco_Human/Core/HOCOMOCOv11_core_HUMAN_mono_jaspar_format_filtered"
        genome="/home/chris/Desktop/Ressources/Human_genome_files/hg38.fa" # genome in fasta format
        #promotor_bed_file=$1 #CSV file from EpiRegio
        #EPIREGIO_output_Dir="./EPIREGIO/" #user defined output dir
        pvalue="0.05" #pvalue threshold for PASTAA
        #folder="dummy"

        #%% Create Output folders
        command="""
        mkdir """ +location+ """/MERILOP/
        mkdir """ +location+ """/MERILOP/input_genes/
        mkdir """ +location+ """/MERILOP/input_2_rest_decreasing/
        mkdir """ +location+ """/MERILOP/input_2_rest_increasing/
        mkdir """ +location+ """/EPIREGIO/
        mkdir """ +location+ """/PASTAA/
        mkdir """ +location+ """/HOMER/
        mkdir """ +location+ """/HOMER_hyper/
        mkdir """ +location+ """/AME_shuffle/
        mkdir """ +location+ """/AME_bg/
        mkdir """ +location+ """/AME_ranking_fisher/
        mkdir """ +location+ """/AME_ranking_3dmhg/
        mkdir """ +location+ """/AME_ranking_4dmhg/
        mkdir """ +location+ """/AME_ranking_fisher_increasing/
        mkdir """ +location+ """/AME_ranking_3dmhg_increasing/
        mkdir """ +location+ """/AME_ranking_4dmhg_increasing/
        mkdir """ +location+ """/AME_ranking_fisher_decreasing/
        mkdir """ +location+ """/AME_ranking_3dmhg_decreasing/
        mkdir """ +location+ """/AME_ranking_4dmhg_decreasing/
        mkdir """ +location+ """/RAMEN/
        """
        for i in command.split("\n"):
            #subprocess.run(i,shell=True)
            pass


    #%%
    command="""
    bedtools getfasta -name -s -fi """+genome+""" -bed """ +location+ """/input_genes_promoters.bed -fo """ +location+ """/input_genes.fa
    bedtools getfasta -name -s -fi """+genome+""" -bed """ +location+ """/input_2rest_decreasing.bed -fo """ +location+ """/input_2rest_decreasing.fa
    bedtools getfasta -name -s -fi """+genome+""" -bed """ +location+ """/input_2rest_increasing.bed -fo """ +location+ """/input_2rest_increasing.fa
    bedtools getfasta -name -s -fi """+genome+""" -bed """ +location+ """/background_promotors.bed -fo """ +location+ """/background_promotors.fa
    
    """

    for i in command.split("\n"):
        subprocess.run(i,shell=True)

    #%% delete the (+)(-) from fasta files gene names because that trips up tools
    with open(location+"/input_genes.fa") as f:
        newText=f.read().replace('(+)', '').replace('(-)', '')

    with open(location+"/input_genes.fa", "w") as f:
        f.write(newText)

    with open(location+"/input_2rest_decreasing.fa") as f:
        newText=f.read().replace('(+)', '').replace('(-)', '')

    with open(location+"/input_2rest_decreasing.fa", "w") as f:
        f.write(newText)

    with open(location+"/input_2rest_increasing.fa") as f:
        newText=f.read().replace('(+)', '').replace('(-)', '')

    with open(location+"/input_2rest_increasing.fa", "w") as f:
        f.write(newText)

    with open(location+"/background_promotors.fa") as f:
        newText=f.read().replace('(+)', '').replace('(-)', '')

    with open(location+"/background_promotors.fa", "w") as f:
        f.write(newText)


    #%% MEIRLOP input file preprocessing
    #AME_REMs_score_adder writes the input text file for MEIRLOP and AME using the promoter names and rankings

    command = """
            echo "MEIRLOP"
            python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ +location+ """/input_genes.fa """ +location+ """/input_genes_ranking_meirlop.txt """ +location+ """/MERILOP/input_genes/ranked_fasta_meirlop_input_genes
            python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ +location+ """/input_2rest_decreasing.fa """ +location+ """/input_2rest_decreasing_meirlop.txt """ +location+ """/MERILOP/input_2_rest_decreasing/ranked_fasta_meirlop_input_2rest_decreasing
            python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ +location+ """/input_2rest_increasing.fa """ +location+ """/input_2rest_increasing_meirlop.txt """ +location+ """/MERILOP/input_2_rest_increasing/ranked_fasta_meirlop_input_2rest_increasing
            """
    bash_out = []
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    # %% MEIRLOP execution

    command = """
                /home/chris/anaconda3/bin/meirlop --fa """ + location + """/MERILOP/input_genes/ranked_fasta_meirlop_input_genes """ + set_of_jaspar_Motifs + """ """ + location + """/MERILOP/input_genes/  
                /home/chris/anaconda3/bin/meirlop --fa """ + location + """/MERILOP/input_2_rest_decreasing/ranked_fasta_meirlop_input_2rest_decreasing """ + set_of_jaspar_Motifs + """ """ + location + """/MERILOP/input_2_rest_decreasing/
                /home/chris/anaconda3/bin/meirlop --fa """ + location + """/MERILOP/input_2_rest_increasing/ranked_fasta_meirlop_input_2rest_increasing """ + set_of_jaspar_Motifs + """ """ + location + """/MERILOP/input_2_rest_increasing/
                """
    bash_out = []
    for i in command.split("\n"):

        #subprocess.run(i,shell=True)

        pass

    #%%
    #PASTAA execution
    command="""
    #PASTAA
    echo "PASTAA"
    echo "PSCM_to_PSEM"
    /home/chris/Desktop/PASTAA/PSCM_to_PSEM """+set_of_transfac_Motifs+""" > """ +location+ """/PASTAA/Energy_Matrix
    
    echo "TRAP input genes"
    /home/chris/Desktop/PASTAA/TRAP """ +location+ """/PASTAA/Energy_Matrix """ +location+ """/input_genes.fa >""" +location+ """/PASTAA/TRAP_output
    echo "Recover_missing_start"
    python3 /home/chris/Desktop/PASTAA/Recover_missing_start.py """ +location+ """/PASTAA/Energy_Matrix """ +location+ """/PASTAA/TRAP_output
    echo "PASTAA"
    /home/chris/Desktop/PASTAA/PASTAA """ +location+ """/PASTAA/trap_recovery """ +location+ """/input_genes_ranking.txt > """ +location+ """/PASTAA/input_ranking
    
    
    echo "TRAP 2 rest decreasing"
    /home/chris/Desktop/PASTAA/TRAP """ +location+ """/PASTAA/Energy_Matrix """ +location+ """/input_2rest_decreasing.fa >""" +location+ """/PASTAA/TRAP_output
    echo "Recover_missing_start"
    python3 /home/chris/Desktop/PASTAA/Recover_missing_start.py """ +location+ """/PASTAA/Energy_Matrix """ +location+ """/PASTAA/TRAP_output
    echo "PASTAA"
    /home/chris/Desktop/PASTAA/PASTAA """ +location+ """/PASTAA/trap_recovery """ +location+ """/input_2rest_decreasing.txt > """ +location+ """/PASTAA/input_2rest_decreasing
    
    echo "TRAP 2 rest increasing"
    /home/chris/Desktop/PASTAA/TRAP """ +location+ """/PASTAA/Energy_Matrix """ +location+ """/input_2rest_increasing.fa >""" +location+ """/PASTAA/TRAP_output
    echo "Recover_missing_start"
    python3 /home/chris/Desktop/PASTAA/Recover_missing_start.py """ +location+ """/PASTAA/Energy_Matrix """ +location+ """/PASTAA/TRAP_output
    echo "PASTAA"
    /home/chris/Desktop/PASTAA/PASTAA """ +location+ """/PASTAA/trap_recovery """ +location+ """/input_2rest_increasing.txt > """ +location+ """/PASTAA/input_2rest_increasing
    
    """
    bash_out=[]
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    #%% AME execution
    #AME initial genes only
    command="""
    #AME shuffle
    echo "AME ranking"
    echo "AME_REMs_score_adder"
    python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ +location+ """/input_genes.fa """ +location+ """/input_genes_ranking.txt """ +location+ """/AME_shuffle/ranked_fasta
    echo "ame"
    /home/chris/Desktop/meme/bin/ame --evalue-report-threshold 10000 --o """+location+"""/AME_shuffle/ame_out/ --control --shuffle-- """ +location+ """/AME_shuffle/ranked_fasta """+set_of_meme_Motifs

    bash_out=[]
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    #AME background
    command="""
    echo "AME ranking"
    echo "AME_REMs_score_adder"
    python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ +location+ """/input_genes.fa """ +location+ """/input_genes_ranking.txt """ +location+ """/AME_bg/ranked_fasta
    echo "ame"
    /home/chris/Desktop/meme/bin/ame --evalue-report-threshold 10000 --o """+location+"""/AME_bg/ame_out/ --control """+location+"""/background_promotors.fa """ +location+ """/AME_bg/ranked_fasta """ +set_of_meme_Motifs

    bash_out=[]
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    #%%
    #AME scored input
    command="""
    echo "AME ranking"
    echo "AME_REMs_score_adder"
    python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """+location+"""/input_genes.fa """+location+"""/input_genes_ranking.txt """+location+"""/AME_ranking_fisher/ranked_fasta
    echo "ame"
    /home/chris/Desktop/meme/bin/ame --evalue-report-threshold 10000 --o """+location+"""/AME_ranking_fisher/ame_out/ """+location+"""/AME_ranking_fisher/ranked_fasta """ +set_of_meme_Motifs+"""
    
    echo "AME decreasing"
    echo "AME_REMs_score_adder2"
    python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_2rest_decreasing.fa """ + location + """/input_2rest_decreasing.txt """ + location + """/AME_ranking_fisher_decreasing/input_2rest_decreasing_ranked_fasta
    echo "AME_REMs_score_adder3"
    echo "ame"
    /home/chris/Desktop/meme/bin/ame --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_fisher_decreasing/ame_out/ """ + location + """/AME_ranking_fisher_decreasing/input_2rest_decreasing_ranked_fasta """ + set_of_meme_Motifs+"""
        
    echo "AME increasing"
    echo "AME_REMs_score_adder"
    python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_2rest_increasing.fa """ + location + """/input_2rest_increasing.txt """ + location + """/AME_ranking_fisher_increasing/input_2rest_increasing_ranked_fasta
    echo "ame"
    /home/chris/Desktop/meme/bin/ame --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_fisher_increasing/ame_out/ """ + location + """/AME_ranking_fisher_increasing/input_2rest_increasing_ranked_fasta """ + set_of_meme_Motifs+"""
        
    """
    bash_out=[]
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    # %%
    # AME 3dmhg
    command = """
            echo "AME ranking"
            echo "AME_REMs_score_adder"
            python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_genes.fa """ + location + """/input_genes_ranking.txt """ + location + """/AME_ranking_3dmhg/ranked_fasta
            echo "ame"
            /home/chris/Desktop/meme/bin/ame --method 3dmhg --scoring totalhits --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_3dmhg/ame_out/ """ + location + """/AME_ranking_3dmhg/ranked_fasta """ + set_of_meme_Motifs+"""
        
            echo "AME ranking"
            echo "AME_REMs_score_adder"
            python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_2rest_decreasing.fa """ + location + """/input_2rest_decreasing.txt """ + location + """/AME_ranking_3dmhg_decreasing/input_2rest_decreasing_ranked_fasta
            echo "ame"
            /home/chris/Desktop/meme/bin/ame --method 3dmhg --scoring totalhits --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_3dmhg_decreasing/ame_out/ """ + location + """/AME_ranking_3dmhg_decreasing/input_2rest_decreasing_ranked_fasta """ + set_of_meme_Motifs+"""
        
            echo "AME ranking"
            echo "AME_REMs_score_adder"
            python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_2rest_increasing.fa """ + location + """/input_2rest_increasing.txt """ + location + """/AME_ranking_3dmhg_increasing/input_2rest_increasing_ranked_fasta
            echo "ame"
            /home/chris/Desktop/meme/bin/ame --method 3dmhg --scoring totalhits --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_3dmhg_increasing/ame_out/ """ + location + """/AME_ranking_3dmhg_increasing/input_2rest_increasing_ranked_fasta """ + set_of_meme_Motifs+"""
        
            """

    bash_out = []
    for i in command.split("\n"):
        #subprocess.run(i, shell=True)
        pass
    # %%
    # AME 4dmhg
    command = """
                echo "AME ranking"
                echo "AME_REMs_score_adder"
                python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_genes.fa """ + location + """/input_genes_ranking.txt """ + location + """/AME_ranking_4dmhg/ranked_fasta
                echo "ame"
                /home/chris/Desktop/meme/bin/ame --method 4dmhg --scoring totalhits --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_4dmhg/ame_out/ """ + location + """/AME_ranking_4dmhg/ranked_fasta """ + set_of_meme_Motifs+"""
                
                echo "AME ranking"
                echo "AME_REMs_score_adder"
                python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_2rest_decreasing.fa """ + location + """/input_2rest_decreasing.txt """ + location + """/AME_ranking_4dmhg_decreasing/input_2rest_decreasing_ranked_fasta
                echo "ame"
                /home/chris/Desktop/meme/bin/ame --method 4dmhg --scoring totalhits --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_4dmhg_decreasing/ame_out/ """ + location + """/AME_ranking_4dmhg_decreasing/input_2rest_decreasing_ranked_fasta """ + set_of_meme_Motifs+"""
        
                echo "AME ranking"
                echo "AME_REMs_score_adder"
                python3 /home/chris/Desktop/meme/AME_REMs_score_adder.py """ + location + """/input_2rest_increasing.fa """ + location + """/input_2rest_increasing.txt """ + location + """/AME_ranking_4dmhg_increasing/input_2rest_increasing_ranked_fasta
                echo "ame"
                /home/chris/Desktop/meme/bin/ame --method 4dmhg --scoring totalhits --evalue-report-threshold 10000 --o """ + location + """/AME_ranking_4dmhg_increasing/ame_out/ """ + location + """/AME_ranking_4dmhg_increasing/input_2rest_increasing_ranked_fasta """ + set_of_meme_Motifs+"""
        
                """

    bash_out = []
    for i in command.split("\n"):
        subprocess.run(i, shell=True)
        #pass
    #%%
    #Ramen execution
    command="""
    echo "RAMEN"
    
    /home/chris/Desktop/meme/libexec/meme-5.1.1/ramen --pvalue-cutoff 1 --bgformat 0  """+location+"""/AME_ranking_fisher/ranked_fasta """+set_of_meme_Motifs+""" >"""+location+"""/RAMEN/RAMEN_output_input_genes
    /home/chris/Desktop/meme/libexec/meme-5.1.1/ramen --pvalue-cutoff 1 --bgformat 0  """+location+"""/AME_ranking_fisher_decreasing/ranked_fasta """+set_of_meme_Motifs+""" >"""+location+"""/RAMEN/RAMEN_output_input_genes_decreasing
    /home/chris/Desktop/meme/libexec/meme-5.1.1/ramen --pvalue-cutoff 1 --bgformat 0  """+location+"""/AME_ranking_fisher_increasing/ranked_fasta """+set_of_meme_Motifs+""" >"""+location+"""/RAMEN/RAMEN_output_input_genes_increasing        
    /home/chris/Desktop/meme/libexec/meme-5.1.1/ramen --pvalue-cutoff 1 --bgformat 0  --bgfile /home/chris/Desktop/Ressources/Analysen/From_chip_to_genelist/promotor_backgrounds_all_genes/500_background.fa """+location+"""/AME_shuffle/ranked_fasta """+set_of_meme_Motifs+""" >"""+location+"""/RAMEN/RAMEN_Output_bg_file
    """

    bash_out=[]
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    #%%
    #HOMER execution
    command="""
    
    echo "HOMER"
    
    
    /home/chris/Desktop/HOMER/bin/findMotifs.pl """+location+"""/input_genes.fa fasta """+location+"""/HOMER/0001/nobg -nomotif -mknown """+set_of_homer_Motifs_0001+"""
    
    /home/chris/Desktop/HOMER/bin/findMotifs.pl """+location+"""/input_genes.fa fasta """+location+"""/HOMER/0001/bg -nomotif -mknown """+set_of_homer_Motifs_0001+""" -fastaBg /home/chris/Desktop/Ressources/Analysen/From_chip_to_genelist/promotor_backgrounds_all_genes/500_background.fa
    
    /home/chris/Desktop/HOMER/bin/findMotifs.pl """+location+"""/input_genes.fa fasta """+location+"""/HOMER/0001/binomial -nomotif -mknown """+set_of_homer_Motifs_0001+""" -b -fastaBg /home/chris/Desktop/Ressources/Analysen/From_chip_to_genelist/promotor_backgrounds_all_genes/500_background.fa
    
    """
    bash_out=[]
    for i in command.split("\n"):
        #subprocess.run(i,shell=True)
        pass

    # %%
    # HOMER execution
    command = """
    
        echo "HOMER"
        
        /home/chris/Desktop/HOMER/bin/findMotifs.pl """ + location + """/input_genes.fa fasta """ + location + """/HOMER/0001/decreasing -nomotif -mknown """ + set_of_homer_Motifs_0001 + """ -b -fastaBg """+location+"""/input_2rest_decreasing.fa 
        /home/chris/Desktop/HOMER/bin/findMotifs.pl """ + location + """/input_genes.fa fasta """ + location + """/HOMER/0001/increasing -nomotif -mknown """ + set_of_homer_Motifs_0001 + """ -b -fastaBg """+location+"""/input_2rest_increasing.fa """

    bash_out = []
    for i in command.split("\n"):
        # subprocess.run(i,shell=True)
        pass

    #%%
    pickle.dump(initial_genesets_going_into_analysis, open("/home/chris/Desktop/Ressources/Analysen/From_chip_to_genelist/Analyse_all_ensemble_promotoren/my_initial_genesets_going_into_analysis_5000.p", "wb"))


