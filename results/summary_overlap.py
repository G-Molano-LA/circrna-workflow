import pandas as pd

filename = "summary_overlap.bed"
counts  = pd.read_csv(filename, sep='\t')

cols = counts.columns # object of type Index that contains column names.
colnames = []
for column in cols:
    colnames.append(column)

colnames_sub = []
colnames_sub = colnames[6:len(colnames)]

count = -1
for name in colnames_sub:
    count = count + 1
    colnames_sub[count] = name.lstrip('/home/alejandra/circrna-workflow/results/identification/overlap/')
    colnames_sub[count] = colnames_sub[count].rstrip('_common.txt_BSJ')

counts = counts.drop(['#chrom','start','end','name','score','strand'], axis = 1)
counts.columns = ['SRR12794681', 'SRR12794682', 'SRR12794683', 'SRR12794684', 'SRR12794685', 'SRR12794686', 'SRR12794687', 'SRR12794688', 'SRR12794689', 'SRR12794690', 'SRR12794691', 'SRR12794692', 'SRR12794693', 'SRR12794694', 'SRR12794695', 'SRR12794696', 'SRR12794697', 'SRR12794698', 'SRR12794699', 'SRR12794700', 'SRR12794701', 'SRR12794702', 'SRR12794703', 'SRR12794704', 'SRR12794705', 'SRR12794706', 'SRR12794707', 'SRR12794708', 'SRR12794709', 'SRR12794710', 'SRR12794711', 'SRR12794712', 'SRR12794713', 'SRR12794714', 'SRR12794715', 'SRR12794716', 'SRR12794717', 'SRR12794718', 'SRR12794719', 'SRR12794720', 'SRR12794721', 'SRR12794722', 'SRR12794723', 'SRR12794724', 'SRR12794725', 'SRR12794726', 'SRR12794727', 'SRR12794728', 'SRR12794729', 'SRR12794730']
