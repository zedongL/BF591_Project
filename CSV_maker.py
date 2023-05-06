import GEOparse
import pandas as pd


# Set working directory
path = r"D:\BF591"
# set file name
file_name = r"\GSE64810_family.soft"

# set file path
soft_file = path + file_name


# load file
gse = GEOparse.get_GEO(filepath=soft_file)
# initialize dataframe
df = pd.DataFrame()
# iterate through metadata
for key, value in gse.gsms.items():
    # iterate through sample metadata
    for k2, v2 in value.metadata.items():
        # if there are multiple values, split them
        if len(v2) > 1:
            for x in range(len(v2)):
                # if there are multiple pairs, split them
                if ":" in v2[x]:
                    k3, v3 = v2[x].split(": ")
                    df.loc[key, k3] = v3
                # if there are multiple values, add them to the same cell
                else:
                    df.loc[key, k2] = " ".join(v2)
        # if there is only one value, add it to the dataframe
        else:
            df.loc[key, k2] = v2[0]

# save dataframe as csv in working directory
df.to_csv(path + file_name.split(".")[0] + ".csv")

print("Job done!")
print("File saved as: " + path + file_name.split(".")[0] + ".csv")
print("information about the file: ")
print(df.info())
