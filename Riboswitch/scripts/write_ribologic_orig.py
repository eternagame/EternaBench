import numpy as np
import pandas as pd
import arnie.utils as utils
import arnie.pfunc as pfunc

def get_lig_constraint(row):
    seq= row['sequence']
    lig = row['ligand']
    lig_aptamer = utils.write_constraints(seq,LIG=True, 
        lig1 = aptamers[lig][0],lig2 = aptamers[lig][1])
    if len(lig_aptamer) > len(seq):
            lig_aptamer = utils.write_constraints(seq, LIG=True,
                lig1 = aptamers['%s_rev' % lig][0],lig2 = aptamers['%s_rev' % lig][1])
    return lig_aptamer


if __name__=='__main__':
	df = pd.read_csv('ribologic_SI.txt',delimiter='\t')

	df = df.drop(['secstruct_ligand','secstruct_noligand'],axis=1)
	df = df.loc[df["ligand"] != 'miRNA']

	#df = df.iloc[100:110] #REMOVE FOR FULL

	df['logkd_nolig'] = np.log(np.where(df['switch']=='OFF', df['Kd_ON'], df['Kd_OFF']))
	df['logkd_lig'] = np.log(np.where(df['switch']=='OFF', df['Kd_OFF'], df['Kd_ON']))
	df['logkd_nolig_scaled'] = df['logkd_nolig'] - np.log(2) # log min value in dataset, k^*
	df['logkd_lig_scaled'] = df['logkd_lig'] - np.log(2) # log min value in dataset, k^*
	df['min_kd_val'] = 2
	df['Round'] = 'Ribologic'
	df['Activation Ratio'] = df['activation_ratio']
	df.replace('Theophylline','Theo',inplace=True)
	df.replace('Tryptophan','Trp',inplace=True)

	df = df.reset_index()

	aptamers = {'FMN':[('nAGGAUAU', '(xxxxxx('),('AGAAGGn', ')xxxxx)')],
	            'FMN_rev':[('AGAAGGn', '(xxxxx('),('nAGGAUAU', ')xxxxxx)')],
	            'Theo':[('GAUACCAG','(xxx(((('),('CCCUUGGCAGC',')xxx)))xxx)')],
	            'Theo_rev':[('CCCUUGGCAGC','(xxx(((xxx('),('GAUACCAG',')xxx))))')],
	            'Trp':[('AGGACCGG','((xxx((('),('CCGCCACU',')))xxx))')],
	            'Trp_rev':[('CCGCCACU','(((xxx(('),('AGGACCGG','))xxx)))')]}

	df['MS2_aptamer'] = df.apply( lambda row: utils.write_constraints(row['sequence'],MS2=True), axis=1)
	df['lig_aptamer'] = df.apply( lambda row: get_lig_constraint(row), axis=1)
	df['MS2_lig_aptamer'] = df.apply( lambda row: utils.write_constraints(row['sequence'], MS2=True, LIG=True, 
	                       lig1 = aptamers[row['ligand']][0],lig2 = aptamers[row['ligand']][1]),axis=1)

	df.to_json('ribologic_orig.json')



