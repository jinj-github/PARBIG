from Bio import SeqIO  
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from Bio import BiopythonParserWarning
import warnings
warnings.filterwarnings("ignore", category=BiopythonParserWarning)

import os
import pandas as pd
import numpy as np
import math

import time
import random
import string

import streamlit as st
import io

import primer3
import base64
import copy

st.set_page_config(
    page_title="PARBiG",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="collapsed"
)

st.markdown("""
<style>
button[kind="header"] {
    display: none !important;
}
.block-container {
    padding-top: 1rem !important;
}
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 style='text-align: left; margin-top: 0px; margin-bottom: -10px; color: #244494;'><a href='https://www.parbig.com' target='_self' style='color: #244494; text-decoration: none;'>PARBiG</a></h1>", unsafe_allow_html=True)

st.html("""<hr style='margin-top: 10px; margin-bottom: 0px;'>""")

col1, col2, col3, col4 = st.columns([10, 1, 6, 3],vertical_alignment="center")

with col1:
    st.image("gene_cluster_scheme.png", caption="Schematic Diagram of Gene Cluster Primer Design")

with col3:
    st.markdown("<h3 style='text-align: left; margin: 0;'>PARBiG</h3>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: left; margin-top: 3px;'>PARBiG is a primer design tool for the cloning and reconstruction of large, complex gene clusters. No programming background is required for users, who can acquire primer design results by uploading the relevant files.</p>", unsafe_allow_html=True)
    
st.html("""<hr style='margin-top: 0px; margin-bottom: 0px;'>""")

file_list = ["Example_BGC.gbk","Example_gene.txt","Example_RBS.txt","Example_vector.gbk"]
text_list = ["Please upload BGC file (.gbk)","Please upload gene list (.txt)","Please upload RBS list (.txt)","Please upload vector backbone file (.gbk)"]
uploaded_files = [] 
for cur_file,cur_text in zip(file_list,text_list):
    with open(cur_file, "rb") as f:
        file_bytes = f.read()
    b64 = base64.b64encode(file_bytes).decode()
    download_link = f'<a href="data:application/octet-stream;base64,{b64}" download="{cur_file}" style="text-decoration: none; color: #244494; margin-left: 15px;font-weight:bold;">Example file</a>'
    st.markdown(
        f'<span style="font-size: 20px; font-weight: bold;">{cur_text}</span>'
        f'{download_link}',
        unsafe_allow_html=True)
    key_name = cur_file.split(".")[0]
    uploaded_file = st.file_uploader(label='data',label_visibility="hidden", key=key_name)
    uploaded_files.append(uploaded_file)
    st.html(f"""<style> .st-key-{key_name} {{margin-top: -25px !important;}} </style>""")

st.markdown('<div style="border: none">', unsafe_allow_html=True)
with st.expander("Advanced Settings", expanded=False):
    st.markdown("<h4 style='text-align: center; margin-top: 10px;'>Primer Design Parameters Configuration</h4>", unsafe_allow_html=True)
    col1, col2, col3, col4 = st.columns(4, gap="medium")

    with col1:
        st.markdown("<p style='text-align: center; font-weight: bold;'>Tm Value (℃)</p>", unsafe_allow_html=True)
        min_tm = st.number_input(
            label="Minimum Tm",
            value=59,  
            min_value=40,
            max_value=80,
            key="min_tm",
            label_visibility="visible"
        )
        opt_tm = st.number_input(
            label="Optimal Tm",
            value=60,  
            min_value=40,
            max_value=80,
            key="opt_tm",
            label_visibility="visible"
        )
        max_tm = st.number_input(
            label="Maximum Tm",
            value=61,  
            min_value=40,
            max_value=80,
            key="max_tm",
            label_visibility="visible"
        )

        max_diff_tm = st.number_input(
            label="Maximum Tm Difference",
            value=1,  
            min_value=1,
            max_value=5,
            key="max_diff_tm",
            label_visibility="visible"
        )

        tm_formula = st.selectbox(
                "Tm Calculation Algorithm",
                options=[0, 1],
                format_func=lambda x: {
                    0: "Breslauer model",
                    1: "Santalucia model",
                }[x],
                index=1,
                help="""PRIMER SALT MONOVALENT is 50 mM;   
                PRIMER SALT DIVALENT is 1.5 mM;  
                PRIMER DNTP CONC is 0.4 mM;  
                PRIMER DNA CONC is 800 nM"""
            )

    with col2:
        st.markdown("<p style='text-align: center; font-weight: bold;'>Primer Length (bp)</p>", unsafe_allow_html=True)
        min_size = st.number_input(
            label="Minimum Length",
            value=18,  
            min_value=15,  
            max_value=29,  
            key="min_size",
            label_visibility="visible"
        )
        opt_size = st.number_input(
            label="Optimal Length",
            value=20,  
            min_value=15,  
            max_value=29,  
            key="opt_size",
            label_visibility="visible"
        )
        max_size = st.number_input(
            label="Maximum Length",
            value=29,  
            min_value=15,  
            max_value=29,  
            key="max_size",
            label_visibility="visible"
        )

    with col3:
        st.markdown("<p style='text-align: center; font-weight: bold;'>GC Content (%)</p>", unsafe_allow_html=True)
        min_gc = st.number_input(
            label="Minimum GC Content",
            value=0,  
            min_value=0,
            max_value=45,
            key="min_gc",
            label_visibility="visible"
        )
        opt_gc = st.number_input(
            label="Optimal GC Content",
            value=50,  
            min_value=min_gc,
            max_value=55,
            key="opt_gc",
            label_visibility="visible"
        )
        max_gc = st.number_input(
            label="Maximum GC Content",
            value=100,  
            min_value=opt_gc,
            max_value=100,
            key="max_gc",
            label_visibility="visible"
        )

    with col4:
        st.markdown("<p style='text-align: center; font-weight: bold;'>Other Parameters</p>", unsafe_allow_html=True)

        vec_insert_pos = st.number_input(
            label="Insertion Site of Gene Cluster in Vector Backbone",
            value=0,  
            key="vec_insert_pos",
            label_visibility="visible"
        )
        homology_arm_len_gv = st.number_input(
            label="Gene-Vector Backbone Homology Arm",
            value=30,  
            key="homology_arm_length_gv",
            label_visibility="visible"
        )
        homology_arm_len_gg = st.number_input(
            label=" Intergenic Homology Arm",
            value=10,  
            key="homology_arm_length_gg",
            label_visibility="visible"
        )
        homology_arm_len_ff = st.number_input(
            label="Intragenic Homology Arm",
            value=40,  
            key="homology_arm_length_ff",
            label_visibility="visible"
        )

        from PIL import Image
        img = Image.open("homology_arm.png") 
        with st.popover("Schematic Diagram"):
            st.image(img, caption="Schematic Diagram", width='stretch')


st.markdown('</div>', unsafe_allow_html=True)

col1, col2, col3 = st.columns([2, 2, 1])
with col2:
    st.markdown("""
                <style>
                div.stButton > button:first-child {
                    background-color: #244494;
                    color:white;
                }
                div.stButton > button:hover {
                    background-color: black;
                    color:white;
                    }
                </style>""", unsafe_allow_html=True)
    run_btn = st.button("Start Primer Design Analysis", key="run_analysis_btn")

st.html("""<hr style='margin-top: 10px; margin-bottom: 0px;'>""")

uploaded_bgc    = uploaded_files[0] #BGC
uploaded_gene   = uploaded_files[1] #gene
uploaded_rbs    = uploaded_files[2] #RBS
uploaded_vector = uploaded_files[3] #vector

if uploaded_bgc is not None and uploaded_gene is not None and uploaded_rbs is not None and uploaded_vector is not None and run_btn:
    # st.write(uploaded_bgc.name,uploaded_gene.name,uploaded_rbs.name,uploaded_vector.name)
    # NRPS11_MG64_BGC57.1a.gbk NRPS11_MG64_BGC57.1a_gene.txt RBS_Strep.txt pSETSc_sfp.gb

    gene = pd.read_table(uploaded_gene,sep="\t")
    gene["Gene"] = gene["Gene"].apply(lambda x: x.strip())

    bgc_name = uploaded_bgc.name.replace(".gbk","")

    genome_list = [ record for record in SeqIO.parse(io.StringIO(uploaded_bgc.read().decode("utf-8")), 'genbank')]
    df = pd.DataFrame()
    for res in genome_list:
        cur = pd.DataFrame([ (
                        feature.qualifiers.get("locus_tag",[None])[0],
                        int(str(feature.location).split(":")[0].replace("[","")),
                        int(str(feature.location).split(":")[1].split("]")[0]),
                        str(feature.location).split("(")[1].replace(")","")
                        )  for feature in res.features if feature.type in ["CDS","gene"] ], columns=["locus_tag", 
                                                                                                        "start", 
                                                                                                        "end", 
                                                                                                        "strand"])
        cur["genome_id"]  = res.id
        cur["genome_seq"] = str(res.seq)

        df = pd.concat([df,cur],ignore_index=True)

    df = df.drop_duplicates()

    def get_gene_seq(x):
        tmp_df = df.loc[df["locus_tag"] == x,:].copy()
        tmp_df["gene_seq"] = tmp_df.apply(lambda m: (m["genome_seq"][ m["start"]:m["end"] ]).translate(str.maketrans("ATCG","TAGC"))[::-1] if m["strand"] == "-" else m["genome_seq"][ m["start"]:m["end"] ], axis=1)
        return tmp_df["gene_seq"].iloc[0]
        
    gene["gene_seq"]        = gene["Gene"].apply(get_gene_seq)
    gene["gene_seq"]        = gene["gene_seq"].apply(lambda x: "ATG" + x[3:])
    gene["gene_seq_length"] = gene["gene_seq"].apply(lambda x: len(x))

    st.markdown('''<h5>Gene sequence</h5>''', unsafe_allow_html=True)
    cur_small = gene.copy().reset_index(drop=True)
    del cur_small["RBS_assignment"]
    st.dataframe(cur_small)

    pre1 = "".join(random.choices(string.ascii_letters, k=10))
    pre2 = time.strftime("%Yy%mm%dd%Hh%Mm%Ss")
    fasta_file = pre1+"_"+pre2+"_"+"generate_gene_seq.fasta"

    with open(fasta_file, "w") as ff:
        for name,seq in zip(gene["Gene"].tolist(),gene["gene_seq"].tolist()):
            ff.write(f'>{name}\n{seq}\n')

    ws = 15
    out = f'{pre1}_{pre2}_blast_{ws}_result.txt'
    command = f"blastn -query {fasta_file} -subject {fasta_file} -dust no -evalue 0.01 -word_size {ws} -outfmt 7 -out {out}"

    os.system(command)

    names = 'query,subject,identity,alignment_length,mismatches,gap_opens,qstart,qend,sstart,send,evalue,bit_score'
    blast_df = pd.read_table(out, sep="\t",comment="#",header=None,names=names.split(","))
    blast_df["evalue"] = blast_df["evalue"].apply(lambda x: round(x,4))
    blast_df["qstart"] = blast_df["qstart"].apply(lambda x: x-1)
    blast_df["sstart"] = blast_df["sstart"].apply(lambda x: x-1)
    blast_df["query_seq_length"] = blast_df["query"].apply(lambda x: gene.loc[gene["Gene"] == x, "gene_seq_length"].iloc[0])

    blast_df = blast_df.loc[~(
                                (blast_df["query"] == blast_df["subject"]) 
                                & (blast_df["identity"] == 100) 
                                & (blast_df["alignment_length"] == blast_df["query_seq_length"])
                            ), :]

    blast75 = blast_df.loc[(blast_df["alignment_length"] > 50) & (blast_df["identity"] > 75),:].sort_values(["identity","alignment_length"],ascending=[False,False])
    blast80 = blast_df.loc[(blast_df["alignment_length"] > 50) & (blast_df["identity"] > 80),:].sort_values(["identity","alignment_length"],ascending=[False,False])
    blast85 = blast_df.loc[(blast_df["alignment_length"] > 50) & (blast_df["identity"] > 85),:].sort_values(["identity","alignment_length"],ascending=[False,False])
    blast90 = blast_df.loc[(blast_df["alignment_length"] > 50) & (blast_df["identity"] > 90),:].sort_values(["identity","alignment_length"],ascending=[False,False])
    
    blast75["repeat_position"] = blast75.apply(lambda x: (x["qstart"], x["qend"]),axis=1)
    group75                    = blast75.groupby("query").agg({"repeat_position":list}).reset_index()
    group75["repeat_position"] = group75["repeat_position"].apply(lambda x: sorted(x))
    
    blast80["repeat_position"] = blast80.apply(lambda x: (x["qstart"], x["qend"]),axis=1)
    group80                    = blast80.groupby("query").agg({"repeat_position":list}).reset_index()
    group80["repeat_position"] = group80["repeat_position"].apply(lambda x: sorted(x))
    
    blast85["repeat_position"] = blast85.apply(lambda x: (x["qstart"], x["qend"]),axis=1)
    group85                    = blast85.groupby("query").agg({"repeat_position":list}).reset_index()
    group85["repeat_position"] = group85["repeat_position"].apply(lambda x: sorted(x))
    
    blast90["repeat_position"] = blast90.apply(lambda x: (x["qstart"], x["qend"]),axis=1)
    group90                    = blast90.groupby("query").agg({"repeat_position":list}).reset_index()
    group90["repeat_position"] = group90["repeat_position"].apply(lambda x: sorted(x))
    
    def merge_overlapping_extended_repeat_position(x):
        if str(x) != "nan":
            sorted_intervals = sorted(x)
            merged = [sorted_intervals[0]]
            for curr_s, curr_e in sorted_intervals[1:]:
                last_s, last_e = merged[-1]
                if curr_s <= last_e:
                    merged[-1] = (last_s, max(last_e, curr_e))
                else:
                    merged.append((curr_s, curr_e))
            return merged
        else:
            return np.nan                                                                               
    
    gene75 = pd.merge(gene,group75,left_on="Gene",right_on="query",how="left")
    del gene75["query"]
    gene75["extended_repeat_position"] = gene75.apply(lambda x: [( max(0,i[0]-400),min(i[1]+400,x["gene_seq_length"]) ) for i in x["repeat_position"]] if str(x["repeat_position"]) != "nan" and x["gene_seq_length"] >= 10000 else np.nan, axis=1)
    gene75["merge_overlapping_extended_repeat_position"] = gene75["extended_repeat_position"].apply(merge_overlapping_extended_repeat_position)
    
    # st.write("=====================75============================")
    # st.write(blast75)
    # st.write(group75)
    # st.write(gene75)
    
    gene80 = pd.merge(gene,group80,left_on="Gene",right_on="query",how="left")
    del gene80["query"]
    gene80["extended_repeat_position"] = gene80.apply(lambda x: [( max(0,i[0]-400),min(i[1]+400,x["gene_seq_length"]) ) for i in x["repeat_position"]] if str(x["repeat_position"]) != "nan" and x["gene_seq_length"] >= 10000 else np.nan, axis=1)
    gene80["merge_overlapping_extended_repeat_position"] = gene80["extended_repeat_position"].apply(merge_overlapping_extended_repeat_position)

    # st.write("=====================80============================")
    # st.write(blast80)
    # st.write(group80)
    # st.write(gene80)
    
    gene85 = pd.merge(gene,group85,left_on="Gene",right_on="query",how="left")
    del gene85["query"]
    gene85["extended_repeat_position"] = gene85.apply(lambda x: [( max(0,i[0]-400),min(i[1]+400,x["gene_seq_length"]) ) for i in x["repeat_position"]] if str(x["repeat_position"]) != "nan" and x["gene_seq_length"] >= 10000 else np.nan, axis=1)
    gene85["merge_overlapping_extended_repeat_position"] = gene85["extended_repeat_position"].apply(merge_overlapping_extended_repeat_position)

    # st.write("=====================85============================")
    # st.write(blast85)
    # st.write(group85)
    # st.write(gene85)
    
    gene90 = pd.merge(gene,group90,left_on="Gene",right_on="query",how="left")
    del gene90["query"]
    gene90["extended_repeat_position"] = gene90.apply(lambda x: [( max(0,i[0]-400),min(i[1]+400,x["gene_seq_length"]) ) for i in x["repeat_position"]] if str(x["repeat_position"]) != "nan" and x["gene_seq_length"] >= 10000 else np.nan, axis=1)
    gene90["merge_overlapping_extended_repeat_position"] = gene90["extended_repeat_position"].apply(merge_overlapping_extended_repeat_position)

    # st.write("=====================90============================")
    # st.write(blast90)
    # st.write(group90)
    # st.write(gene90)
    
    judge75 = gene75["merge_overlapping_extended_repeat_position"].apply(lambda x: any(i[1]-i[0] >= 10000 for i in x) if str(x) != "nan" else False).any()
    judge80 = gene80["merge_overlapping_extended_repeat_position"].apply(lambda x: any(i[1]-i[0] >= 10000 for i in x) if str(x) != "nan" else False).any()
    judge85 = gene85["merge_overlapping_extended_repeat_position"].apply(lambda x: any(i[1]-i[0] >= 10000 for i in x) if str(x) != "nan" else False).any()
    judge90 = gene90["merge_overlapping_extended_repeat_position"].apply(lambda x: any(i[1]-i[0] >= 10000 for i in x) if str(x) != "nan" else False).any()

    # st.write(judge75)    
    # st.write(judge80)    
    # st.write(judge85)    
    # st.write(judge90)    
    
    if not judge75:
        del gene
        gene = gene75.copy()
        with st.expander("Sequence alignment", expanded=False):
            cur_small = blast75.copy().reset_index(drop=True)
            if cur_small.shape[0] != 0:
                st.dataframe(cur_small)
            else:
                st.markdown('''<h6>No significant repetitive sequences</h6>''', unsafe_allow_html=True)
    elif not judge80:
        del gene
        gene = gene80.copy()
        with st.expander("Sequence alignment", expanded=False):
            cur_small = blast80.copy().reset_index(drop=True)
            if cur_small.shape[0] != 0:
                st.dataframe(cur_small)
            else:
                st.markdown('''<h6>No significant repetitive sequences</h6>''', unsafe_allow_html=True)
    elif not judge85:
        del gene
        gene = gene85.copy()
        with st.expander("Sequence alignment", expanded=False):
            cur_small = blast85.copy().reset_index(drop=True)
            if cur_small.shape[0] != 0:
                st.dataframe(cur_small)
            else:
                st.markdown('''<h6>No significant repetitive sequences</h6>''', unsafe_allow_html=True)
    elif not judge90:
        del gene
        gene = gene90.copy()
        with st.expander("Sequence alignment", expanded=False):
            cur_small = blast90.copy().reset_index(drop=True)
            if cur_small.shape[0] != 0:
                st.dataframe(cur_small)
            else:
                st.markdown('''<h6>No significant repetitive sequences</h6>''', unsafe_allow_html=True)
    else:
        st.warning("The length of the repeat region is ≥ 10000!")
        st.stop()

    def get_max_segment_length(x):
        if x["gene_seq_length"] <= 10000:
            return np.nan
        else:
            tmp1 = math.ceil(x["gene_seq_length"] / 10000)
            tmp2 = min(9000,math.floor(x["gene_seq_length"] / tmp1))
            return tmp2

    gene["max_segment_length"] = gene.apply(get_max_segment_length,axis=1)   

    def get_non_repeat_position(x):
        if str(x["merge_overlapping_extended_repeat_position"]) != "nan" and str(x["max_segment_length"]) != "nan":
            non_rep = []
            last = 0
            for s, e in x["merge_overlapping_extended_repeat_position"]:
                if s > last:
                    non_rep.append((last, s))
                last = e
            if last < x["gene_seq_length"]:
                non_rep.append((last, x["gene_seq_length"]))
            return non_rep
        else:
            if pd.isna(x["merge_overlapping_extended_repeat_position"]) and str(x["max_segment_length"]) != "nan":
                num = math.ceil(x["gene_seq_length"] / x["max_segment_length"])
                seq_list = [(i * x["max_segment_length"], min((i + 1) * x["max_segment_length"], x["gene_seq_length"])) for i in range(num)]
                return seq_list
            else:
                return np.nan
                
    gene["non_repeat_position"] = gene.apply(get_non_repeat_position,axis=1)

    def get_small_units(x):
        if str(x["non_repeat_position"]) != "nan" and str(x["merge_overlapping_extended_repeat_position"]) != "nan":
            units = []
            for s, e in x["non_repeat_position"]:
                pos = s
                while pos < e:
                    end = min(pos + 10, e)
                    units.append((pos, end))
                    pos = end
            for s, e in x["merge_overlapping_extended_repeat_position"]:
                units.append((s, e))
            units = sorted(units)
            return units
        else:
            if str(x["non_repeat_position"]) != "nan" and pd.isna(x["merge_overlapping_extended_repeat_position"]):
                return x["non_repeat_position"]
            else:
                return np.nan

    gene["small_units"] = gene.apply(get_small_units,axis=1)

    def merge_small_units(x):
        if str(x["small_units"]) != "nan" and str(x["max_segment_length"]) != "nan":
            
            units = x["small_units"]
            max_len = x["max_segment_length"]
            segs = []
            curr_start, cum_len = units[0][0], 0
        
            for s, e in units:
                l = e - s  
                if cum_len + l > max_len:
                    segs.append((curr_start, s))
                    curr_start = s
                    cum_len = l
                else:
                    cum_len += l
        
            segs.append((curr_start, units[-1][1]))

            while len(segs) >= 2:
                idx = len(segs) - 1
                cur_len = segs[idx][1] - segs[idx][0]
                if cur_len < 1000:
                    segs[idx-1] = (segs[idx-1][0], segs[idx][1])
                    segs.pop(idx)
                else:
                    break

            return segs
            
        else:
            return np.nan

    gene["merge_small_units"] = gene.apply(merge_small_units,axis=1)

    def processed_segment(x):
        if str(x) != "nan":
            act_num  = len(x)
            theo_num = math.ceil(x[-1][1] / 10000)
            
            if act_num == theo_num:
                return x
                
            diff_num = act_num - theo_num
            segments = copy.deepcopy(x)
            for _ in range(diff_num):
                merge_lengths = sorted([ (segments[i+1][1] - segments[i][0],i) for i in range(len(segments) - 1) ])
                min_idx = merge_lengths[0][1]
                merged = (segments[min_idx][0], segments[min_idx+1][1])
                segments = segments[:min_idx] + [merged] + segments[min_idx+2:]
            return segments
        else:
            return np.nan
    
    gene["merge_small_units"] = gene["merge_small_units"].apply(processed_segment)       
                
    if gene["merge_small_units"].apply(lambda x: str(x) != "nan").any():
        st.markdown('''<h5>Gene sequence segmentation</h5>''', unsafe_allow_html=True)
        cur_small = gene.loc[:,["Gene","gene_seq_length","repeat_position","merge_small_units"]].reset_index(drop=True)
        cur_small.columns = ["Gene","gene_seq_length","repetitive_position","seq_segment_region"]
        cur_small["repetitive_position"] = cur_small["repetitive_position"].astype(str)
        cur_small["seq_segment_region"] = cur_small["seq_segment_region"].astype(str)
        st.dataframe(cur_small)

    # st.write(min_size)
    # st.write(opt_size)
    # st.write(max_size)
    # st.write(min_tm)
    # st.write(opt_tm)
    # st.write(max_tm)
    # st.write(min_gc)
    # st.write(opt_gc)
    # st.write(max_gc)
    # st.write(tm_formula)
    # st.write(max_diff_tm)    

    import primer_design  

    gene_segs_overlap_length = homology_arm_len_ff  
    def get_raw_clone_primer(x):
        if str(x["merge_small_units"]) != "nan":
            cur_pos_list = []
            cur_pos_list.append((x["merge_small_units"][0][0], x["merge_small_units"][0][1]+gene_segs_overlap_length))
            for middle_pos in x["merge_small_units"][1:-1]:
                cur_pos_list.append((middle_pos[0]-gene_segs_overlap_length, middle_pos[1]+gene_segs_overlap_length))
            cur_pos_list.append((x["merge_small_units"][-1][0]-gene_segs_overlap_length, x["merge_small_units"][-1][1]))
            cur_seq_list = [x["gene_seq"][ int(i[0]):int(i[1]) ] for i in cur_pos_list]
            cur_primer_list = [primer_design.get_primer3_data_for_clone(i,min_size,opt_size,max_size,min_tm,opt_tm,max_tm,min_gc,opt_gc,max_gc,tm_formula,max_diff_tm) for i in cur_seq_list]
            return cur_primer_list
        else:
            return [primer_design.get_primer3_data_for_clone(x["gene_seq"],min_size,opt_size,max_size,min_tm,opt_tm,max_tm,min_gc,opt_gc,max_gc,tm_formula,max_diff_tm)]
        
    gene["raw_clone_primer"] = gene.apply(get_raw_clone_primer,axis=1)

    # st.write(vec_insert_pos)
    # st.write(homology_arm_len_gv)
    # st.write(homology_arm_len_gg)
    # st.write(homology_arm_len_ff)

    vector = [ record for record in SeqIO.parse(io.StringIO(uploaded_vector.read().decode("utf-8")), 'genbank')][0]
    vector_seq = str(vector.seq).upper()  

    # st.write(vector_seq[0:10])
    # st.write(vector_seq[-10:])

    vector_insert_pos = vec_insert_pos
    gene_vector_overlap_length = homology_arm_len_gv 

    if vector_insert_pos < gene_vector_overlap_length:
        dif = gene_vector_overlap_length - vector_insert_pos
        first_gene_vector_overlap_seq = vector_seq[(-1)*dif:] + vector_seq[0:vector_insert_pos]
        last_gene_vector_overlap_seq  = vector_seq[vector_insert_pos:vector_insert_pos+gene_vector_overlap_length]
    else:
        first_gene_vector_overlap_seq = vector_seq[vector_insert_pos - gene_vector_overlap_length:vector_insert_pos]
        if len(vector_seq) - vector_insert_pos > gene_vector_overlap_length:
            last_gene_vector_overlap_seq  = vector_seq[vector_insert_pos:vector_insert_pos+gene_vector_overlap_length]
        else:
            dif = gene_vector_overlap_length - (len(vector_seq) - vector_insert_pos)
            last_gene_vector_overlap_seq  = vector_seq[vector_insert_pos:] +  vector_seq[0:dif]        
    
    gene_number = gene.shape[0]
    gene["up_vector_overlap"]   =  [first_gene_vector_overlap_seq] + [ np.nan for i in range(gene_number - 1) ]
    gene["down_vector_overlap"] =  [ np.nan for i in range(gene_number - 1) ] + [last_gene_vector_overlap_seq]   

    gene_rbs_overlap_length = homology_arm_len_gg 

    if gene.shape[0] != 1:
        gene["up_gene_overlap"]   =  [np.nan] + [i[(-1)*gene_rbs_overlap_length:] for i in gene["gene_seq"].tolist()[0:-1]]
        gene["down_gene_overlap"] =  [i[0:gene_rbs_overlap_length] for i in gene["gene_seq"].tolist()[1:]] + [np.nan]    
    else:
        gene["up_gene_overlap"]   = [np.nan]
        gene["down_gene_overlap"] = [np.nan]

    if gene.shape[0] != 1:
        rbs = pd.read_table(uploaded_rbs,sep="\t")
        rbs["RBS"] = rbs["RBS"].apply(lambda x: x.strip().upper())
        gene["up_RBS"]  = [np.nan] + rbs["RBS"].tolist()[0:gene_number - 1]
        gene["down_RBS"] = rbs["RBS"].tolist()[0:gene_number - 1] + [np.nan]    
    else:
        gene["up_RBS"]   = [np.nan]
        gene["down_RBS"] = [np.nan]
        
    gene["up_vector_overlap"]   =  gene["up_vector_overlap"].fillna('')
    gene["down_vector_overlap"] =  gene["down_vector_overlap"].fillna('')
    gene["up_gene_overlap"]     =  gene["up_gene_overlap"].fillna('')
    gene["down_gene_overlap"]   =  gene["down_gene_overlap"].fillna('')
    gene["up_RBS"]              =  gene["up_RBS"].fillna('')
    gene["down_RBS"]            =  gene["down_RBS"].fillna('')    

    def get_processed_clone_primer_seq(x):
        if len(x["raw_clone_primer"]) > 1:
            processed_list = []
            processed_list.append( (x["up_vector_overlap"] + x["up_gene_overlap"] + x["up_RBS"] + x["raw_clone_primer"][0][0], x["raw_clone_primer"][0][1]) )
            for cur in x["raw_clone_primer"][1:-1]:
                processed_list.append( (cur[0],cur[1]) )
            processed_list.append( (x["raw_clone_primer"][-1][0], x["down_vector_overlap"].translate(str.maketrans("ATCG","TAGC"))[::-1] + x["down_gene_overlap"].translate(str.maketrans("ATCG","TAGC"))[::-1] + x["down_RBS"].translate(str.maketrans("ATCG","TAGC"))[::-1] + x["raw_clone_primer"][-1][1]) )
            return processed_list
        else:
            processed_list = []
            processed_list.append( (x["up_vector_overlap"] + x["up_gene_overlap"] + x["up_RBS"] + x["raw_clone_primer"][0][0], x["down_vector_overlap"].translate(str.maketrans("ATCG","TAGC"))[::-1] + x["down_gene_overlap"].translate(str.maketrans("ATCG","TAGC"))[::-1] + x["down_RBS"].translate(str.maketrans("ATCG","TAGC"))[::-1] + x["raw_clone_primer"][0][1]) )
            return processed_list

    gene["processed_clone_primer_seq"] = gene.apply(get_processed_clone_primer_seq,axis=1)   

    gene["raw_processed"] = gene.apply(lambda x:[(m,n) for m,n in zip(x["raw_clone_primer"],x["processed_clone_primer_seq"])],axis=1)

    vector_clone_overlap_gene_length = 10  

    new_vector_seq = vector_seq[vector_insert_pos:] + vector_seq[0:vector_insert_pos]

    raw_vector_primer = primer_design.get_primer3_data_for_clone(new_vector_seq,min_size,opt_size,max_size,min_tm,opt_tm,max_tm,min_gc,opt_gc,max_gc,tm_formula,max_diff_tm)

    processed_vector_primer = ( 
                                gene["gene_seq"].iloc[-1][(-1)*vector_clone_overlap_gene_length:] + raw_vector_primer[0],
                                (gene["gene_seq"].iloc[0][0:vector_clone_overlap_gene_length]).translate(str.maketrans("ATCG","TAGC"))[::-1] + raw_vector_primer[1]
                            )     

    clone_list = [ (raw_vector_primer, processed_vector_primer) ] + [ i for cur_list in gene["raw_processed"].tolist() for i in cur_list]
    clone_df   = pd.DataFrame({"raw_seq": [ i[0] for i in clone_list ],"processed_seq": [ i[1] for i in clone_list ]})   

    clone_df["num"] = [ str(i) for i in range(clone_df.shape[0])]

    clone_df["Forward primer name"]     = clone_df["num"].apply(lambda x: f"{bgc_name}_F{x}F")
    clone_df["Forward primer sequence"] = clone_df["processed_seq"].apply(lambda x: x[0])
    clone_df["Forward primer Tm"]       = clone_df["raw_seq"].apply(lambda x: x[4])

    clone_df["Reverse primer name"]     = clone_df["num"].apply(lambda x: f"{bgc_name}_F{x}R")
    clone_df["Reverse primer sequence"] = clone_df["processed_seq"].apply(lambda x: x[1])
    clone_df["Reverse primer Tm"]       = clone_df["raw_seq"].apply(lambda x: x[5])

    st.markdown('''<h5>Amplification primers</h5>''', unsafe_allow_html=True)
    cur_small = clone_df.iloc[:,3:].reset_index(drop=True)
    st.dataframe(cur_small)

    check_point_length = 400

    def get_intragenic_check_seq(x):
        if str(x["merge_small_units"]) != "nan":
            point_list = [end for start, end in x["merge_small_units"][:-1]]
            return [ x["gene_seq"][int(cur - check_point_length) : int(cur + check_point_length)] for cur in point_list ]
        else:
            return [np.nan]

    gene["intragenic_check_seq"] = gene.apply(get_intragenic_check_seq,axis=1) 
    
    seq_list = gene["gene_seq"].tolist()
    if len(seq_list) == 1:
        gene["intergenic_check_seq"] = np.nan
    else:
        gene["intergenic_check_seq"] = [seq_list[i][(-1)*check_point_length:] + seq_list[i+1][0:check_point_length] for i in range(len(seq_list)-1)] + [np.nan]

    tmp = []
    for m,n in zip(gene["intragenic_check_seq"].tolist(),gene["intergenic_check_seq"]):
        tmp.extend(m)
        tmp.append(n)
    new = [ i for i in tmp if pd.notna(i)]

    first_check_seq = new_vector_seq[(-1)*check_point_length:] + gene["gene_seq"].iloc[0][0:check_point_length]
    last_check_seq  = gene["gene_seq"].iloc[-1][(-1)*check_point_length:] + new_vector_seq[0:check_point_length]

    # st.write(min_size)
    # st.write(opt_size)
    # st.write(max_size)
    # st.write(min_tm)
    # st.write(opt_tm)
    # st.write(max_tm)
    # st.write(min_gc)
    # st.write(opt_gc)
    # st.write(max_gc)
    # st.write(tm_formula)
    # st.write(max_diff_tm)    

    all_check_seq  = [first_check_seq] + new + [last_check_seq]

    check_seq_df = pd.DataFrame({"seq":all_check_seq})
    check_seq_df["primer"] = check_seq_df["seq"].apply(lambda x: primer_design.get_primer3_data_for_detection(x,min_size,opt_size,max_size,min_tm,opt_tm,max_tm,min_gc,opt_gc,max_gc,tm_formula,max_diff_tm))

    check_seq_df["num"] = [ str(i+1) for i in range(check_seq_df.shape[0])]
    check_seq_df["Forward primer name"]     = check_seq_df["num"].apply(lambda x: f"{bgc_name}_J{x}F")
    check_seq_df["Forward primer sequence"] = check_seq_df["primer"].apply(lambda x: x[0])
    check_seq_df["Forward primer Tm"]       = check_seq_df["primer"].apply(lambda x: x[4])

    check_seq_df["Reverse primer name"]     = check_seq_df["num"].apply(lambda x: f"{bgc_name}_J{x}R")
    check_seq_df["Reverse primer sequence"] = check_seq_df["primer"].apply(lambda x: x[1])
    check_seq_df["Reverse primer Tm"]       = check_seq_df["primer"].apply(lambda x: x[5])    
 
    st.markdown('''<h5>Verification primers</h5>''', unsafe_allow_html=True)
    cur_small = check_seq_df.iloc[:,3:].reset_index(drop=True)

    for i in cur_small.columns.tolist():
         cur_small[i] = cur_small[i].astype(str)

    st.dataframe(cur_small)

    if gene.shape[0] != 1:
        seq_list_a = [(f"gene_{x}",y) for x,y in zip(gene["Gene"].tolist(),gene["gene_seq"].tolist())]
        num = len(seq_list_a)
        seq_list_b = [(f"RBS_{x}",y) for x,y in enumerate(rbs["RBS"].tolist()[0:num - 1])]
        merged = [elem for pair in zip(seq_list_a, seq_list_b + [None]) for elem in pair if elem is not None]
    else:
        merged = [(f"gene_{x}",y) for x,y in zip(gene["Gene"].tolist(),gene["gene_seq"].tolist())]

    merged.append(("vector",new_vector_seq))

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    total_seq = ""
    elem_positions = []
    for elem_name, elem_seq in merged:
        start = len(total_seq)
        end = start + len(elem_seq)
        elem_positions.append((elem_name, start, end))
        total_seq += elem_seq

    record = SeqRecord(
        seq=Seq(total_seq),  
        id="Custom_Vector",
        name="Custom",
        description="Simple circular vector",
        annotations={
            "molecule_type": "DNA",  
            "topology": "circular"
        }
    )

    features = []
    for elem_name, start, end in elem_positions:
        if elem_name.startswith("gene_"):
            feat_type = "gene"
            qualifiers = {"label": elem_name}
        elif elem_name.startswith("RBS_"):
            feat_type = "misc_feature"
            qualifiers = {"label": elem_name,"note": "Ribosome Binding Site (RBS)"}
        else:
            feat_type = "misc_feature"
            qualifiers = {"label": elem_name,"note": "Vector Backbone (circular)"}
        
        features.append(SeqFeature(
            FeatureLocation(start, end),
            type=feat_type,
            qualifiers=qualifiers
        ))

    record.features = features

    output_path = bgc_name + "_" + pre2 + "_assembly.gbk"
    SeqIO.write(record, output_path, "genbank")    

    with open(output_path, "r", encoding="utf-8") as f:
            gbk_content = f.read()

    st.download_button(
        label="📥 Download assembly file",
        data=gbk_content,
        file_name=os.path.basename(output_path),  
        mime="text/plain",
        key="download_local_gbk"
    )
