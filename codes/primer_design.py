import primer3
import pandas as pd

def get_primer3_data_for_clone(x,
                               min_size,opt_size,max_size,
                               min_tm,opt_tm,max_tm,
                               min_gc,opt_gc,max_gc,
                               tm_formula,max_diff_tm
                              ):

    seq_args = {
                'SEQUENCE_ID' : "cur",
                'SEQUENCE_TEMPLATE': x,
                'SEQUENCE_INCLUDED_REGION': [0,len(x)],
                'SEQUENCE_FORCE_LEFT_START': 0,
                'SEQUENCE_FORCE_RIGHT_START': len(x)-1
                }
    global_args = {
                'PRIMER_TASK':'generic',
        
                'PRIMER_PICK_LEFT_PRIMER':1,
                'PRIMER_PICK_INTERNAL_OLIGO':0,
                'PRIMER_PICK_RIGHT_PRIMER':1,
        
                'PRIMER_PRODUCT_SIZE_RANGE':[30,len(x)],
        
                'PRIMER_MIN_SIZE':min_size,
                'PRIMER_OPT_SIZE':opt_size,
                'PRIMER_MAX_SIZE':max_size,

                'PRIMER_TM_FORMULA':tm_formula,  
                
                "PRIMER_SALT_MONOVALENT": 50.0,  
                "PRIMER_SALT_DIVALENT": 1.5,     
                "PRIMER_DNTP_CONC": 0.4,         
                "PRIMER_DNA_CONC": 800.0,        
                "PRIMER_SALT_CORRECTIONS": 1,    
                
                'PRIMER_MIN_TM':min_tm,
                'PRIMER_OPT_TM':opt_tm,
                'PRIMER_MAX_TM':max_tm,

                'PRIMER_PAIR_MAX_DIFF_TM':max_diff_tm,

                'PRIMER_MIN_GC':min_gc,
                'PRIMER_OPT_GC_PERCENT':opt_gc,
                'PRIMER_MAX_GC':max_gc,

                'PRIMER_NUM_RETURN':1,

                'PRIMER_WT_TM_LT':10,
                'PRIMER_WT_TM_GT':10,

                "PRIMER_GC_CLAMP":0,
                "PRIMER_MAX_SELF_ANY":100,
                "PRIMER_MAX_SELF_END":100,
                "PRIMER_MAX_HAIRPIN_TH":100,
                "PRIMER_MAX_POLY_X":100,
                "PRIMER_MAX_END_STABILITY":100,
        
                'PRIMER_PICK_ANYWAY': 1
        
                }
    
    primer3_result = primer3.bindings.design_primers(seq_args, global_args)

    res = (primer3_result.get('PRIMER_LEFT_0_SEQUENCE',None),
            primer3_result.get('PRIMER_RIGHT_0_SEQUENCE',None),
            primer3_result.get('PRIMER_LEFT_0',None),
            primer3_result.get('PRIMER_RIGHT_0',None),
            primer3_result.get('PRIMER_LEFT_0_TM',None),
            primer3_result.get('PRIMER_RIGHT_0_TM',None),
            primer3_result.get('PRIMER_LEFT_0_GC_PERCENT',None),
            primer3_result.get('PRIMER_RIGHT_0_GC_PERCENT',None))

    return res  


def get_primer3_data_for_detection(x,
                               min_size,opt_size,max_size,
                               min_tm,opt_tm,max_tm,
                               min_gc,opt_gc,max_gc,
                               tm_formula,max_diff_tm
                              ):

    seq_args = {
                'SEQUENCE_ID' : "cur",
                'SEQUENCE_TEMPLATE': x,
                'SEQUENCE_INCLUDED_REGION': [0,len(x)],

                'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0,100,len(x)-100,100]
        
                }
    global_args = {
                'PRIMER_TASK':'generic',
        
                'PRIMER_PICK_LEFT_PRIMER':1,
                'PRIMER_PICK_INTERNAL_OLIGO':0,
                'PRIMER_PICK_RIGHT_PRIMER':1,
        
                'PRIMER_PRODUCT_SIZE_RANGE':[30,len(x)],
        
                'PRIMER_MIN_SIZE':min_size,
                'PRIMER_OPT_SIZE':opt_size,
                'PRIMER_MAX_SIZE':max_size,

                'PRIMER_TM_FORMULA':tm_formula,  
                
                "PRIMER_SALT_MONOVALENT": 50.0,  
                "PRIMER_SALT_DIVALENT": 1.5,     
                "PRIMER_DNTP_CONC": 0.4,         
                "PRIMER_DNA_CONC": 800.0,        
                "PRIMER_SALT_CORRECTIONS": 1,    
                
                'PRIMER_MIN_TM':min_tm,
                'PRIMER_OPT_TM':opt_tm,
                'PRIMER_MAX_TM':max_tm,

                'PRIMER_PAIR_MAX_DIFF_TM':max_diff_tm,

                'PRIMER_MIN_GC':min_gc,
                'PRIMER_OPT_GC_PERCENT':opt_gc,
                'PRIMER_MAX_GC':max_gc,

                'PRIMER_NUM_RETURN':1,
        
                'PRIMER_WT_TM_LT':10,
                'PRIMER_WT_TM_GT':10,

                "PRIMER_GC_CLAMP":0,
                "PRIMER_MAX_SELF_ANY":100,
                "PRIMER_MAX_SELF_END":100,
                "PRIMER_MAX_HAIRPIN_TH":100,
                "PRIMER_MAX_POLY_X":100,
                "PRIMER_MAX_END_STABILITY":100,
        
                'PRIMER_PICK_ANYWAY': 1
        
                }
    
    primer3_result = primer3.bindings.design_primers(seq_args, global_args)

    print(primer3_result)

    res = (primer3_result.get('PRIMER_LEFT_0_SEQUENCE',None),
            primer3_result.get('PRIMER_RIGHT_0_SEQUENCE',None),
            primer3_result.get('PRIMER_LEFT_0',None),
            primer3_result.get('PRIMER_RIGHT_0',None),
            primer3_result.get('PRIMER_LEFT_0_TM',None),
            primer3_result.get('PRIMER_RIGHT_0_TM',None),
            primer3_result.get('PRIMER_LEFT_0_GC_PERCENT',None),
            primer3_result.get('PRIMER_RIGHT_0_GC_PERCENT',None))

    return res  





