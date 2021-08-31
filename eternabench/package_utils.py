
ARNIEDIR = '/home/users/hannahw1/arnie'
EFOLD_PARAMS_LOC = '/home/users/hannahw1/efold_params'

package_options = {
'vienna_2': {'package': 'vienna_2'},
'vienna_2_nodangles': {'package':'vienna_2', 'dangles': False},
'vienna_2_60C': {'package': 'vienna_2', 'T': 60},
'vienna_1': {'package': 'vienna_1'},

'nupack_99': {'package':'nupack'},
'nupack_99_60C': {'package':'nupack','T':60},
'nupack_99_nodangles':{'package':'nupack', 'dangles': False},
'nupack_95':{'package':'nupack_95'},
'nupack_95_nodangles':{'package':'nupack_95', 'dangles': False},

'rnastructure': {'package':'rnastructure'},
'rnastructure_nocoax': {'package':'rnastructure', 'coaxial':False},
'rnastructure_60C': {'package':'rnastructure','T':60},

'contrafold_1': {'package':'contrafold_1'},
'contrafold_2': {'package':'contrafold_2'},
'contrafold_2_nc': {'package':'contrafold_2','param_file':'%s/parameter_files/contrafold.params.noncomplementary' % ARNIEDIR},

'vienna_langdon_pars': {'package': 'vienna_2', 'param_file':'%s/parameter_files/rna_langdon2018.par' % ARNIEDIR},
'vienna_rnasoft_pars': {'package': 'vienna_2', 'param_file':'%s/parameter_files/rna_andronescu2007.par' % ARNIEDIR},
'rnasoft_99' : {'package':'rnasoft_99'},
'rnasoft_07' : {'package':'rnasoft_07'},
'rnasoft_blstar' : {'package':'rnasoft'},
'rnasoft_99_nodangles' : {'package':'rnasoft_99-no-dangles'},
'rnasoft_bl_nodangles' : {'package':'rnasoft_bl-no-dangles'},
'rnasoft_lam-cg' : {'package':'rnasoft_lam-cg'},
'rnasoft_nom-cg' : {'package':'rnasoft_nom-cg'},

'eternafold_A_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_A' % EFOLD_PARAMS_LOC},
'eternafold_B_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_B' % EFOLD_PARAMS_LOC},
'eternafold_C_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_C' % EFOLD_PARAMS_LOC},
'eternafold_D_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_D' % EFOLD_PARAMS_LOC},
'eternafold_E_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_E' % EFOLD_PARAMS_LOC},
'eternafold_F_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_F' % EFOLD_PARAMS_LOC},
'eternafold_G_cfold': {'package': 'contrafold_2', 'param_file':'%s/params_G' % EFOLD_PARAMS_LOC},

'eternafold_A': {'package': 'eternafold', 'param_file':'%s/params_A' % EFOLD_PARAMS_LOC},
'eternafold': {'package': 'eternafold'},
'eternafold_B': {'package': 'eternafold'},
'eternafold_C': {'package': 'eternafold', 'param_file':'%s/params_C' % EFOLD_PARAMS_LOC},
'eternafold_D': {'package': 'eternafold', 'param_file':'%s/params_D' % EFOLD_PARAMS_LOC},
'eternafold_E': {'package': 'eternafold', 'param_file':'%s/params_E' % EFOLD_PARAMS_LOC},
'eternafold_F': {'package': 'eternafold', 'param_file':'%s/params_F' % EFOLD_PARAMS_LOC},
'eternafold_G': {'package': 'eternafold', 'param_file':'%s/params_G' % EFOLD_PARAMS_LOC},

 }

example_subset = ['vienna_2', 'vienna_2_60C', 'nupack_99','rnastructure',
'rnastructure_60C', 'contrafold_2', 'rnasoft_blstar','eternafold_B']

