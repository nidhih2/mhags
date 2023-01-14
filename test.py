import pandas as pd 


x = pd.read_csv("donor_and_host_data_merged.txt", sep="\t")
xx = x[["host_GT", "donor_GT","host_REF","host_ALT", "donor_REF","donor_ALT"]]
xx.to_csv("just5.txt", sep="\t", index=False)

datatemp = xx


datatemp = datatemp[((datatemp['Host_1']!=datatemp['Donor_1']) & (datatemp['Host_1']!=datatemp['Donor_2'])) | 
				((datatemp['Host_2']!=datatemp['Donor_1']) & (datatemp['Host_2']!=datatemp['Donor_2']))]
