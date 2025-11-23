import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def scatter(df,xlabel,ylabel,out,highlight_genes=None,label_offset=0.05):
    indicator = ((df >= 1).sum(axis=1) >= 2)
    df = df[indicator]
    fig,ax = plt.subplots(figsize=(20,20))
    sns.scatterplot(
        x=np.log(1+df.iloc[:,:2].mean(axis=1)),y=np.log(1+df.iloc[:,2:].mean(axis=1)),ax=ax,color='lightgray'
        )
    # 如果有需要高亮的基因
    if highlight_genes is not None:
        highlight_df = df.loc[highlight_genes]
        if not highlight_df.empty:
            x_highlight = np.log(1 + highlight_df.iloc[:, :2].mean(axis=1))
            y_highlight = np.log(1 + highlight_df.iloc[:, 2:].mean(axis=1))
            sns.scatterplot(
                x=x_highlight, y=y_highlight, ax=ax, color='red'
            )
        for gene, x_pos, y_pos in zip(highlight_genes, x_highlight, y_highlight):
            ax.text(x_pos + label_offset, y_pos + label_offset, gene,fontsize=12,color='black')
    # 添加坐标轴标签
    ax.set_xlabel(f"log1p {xlabel} (counts per million)")
    ax.set_ylabel(f"log1p {ylabel} (counts per million)")
    plt.savefig(f"{out}.png",dpi=300)
    plt.close(fig)

# --------------------------
df = pd.read_csv("SLAMseq.new.csv",header=0,index_col=0)
df = df.join(pd.read_csv("Gene.info.csv",index_col=0).iloc[:,0]).set_index("Symbol")
# cpm normalization
cpm = df.apply(lambda x:x*1000000/x.sum(),axis="index")

c1c2 = cpm[["C2.1","C2.2","C1.1","C1.2"]]
c3c4 = cpm[["C4.1","C4.2","C3.1","C3.2"]]

scatter(c1c2,"mEGFP-Rapa-2d","mEGFP-Rapa-2d","SLAMseq.new.Rapa",highlight_genes=["CDK2","HBZ","ESPN","HBG2","HBG1","ENSG00000284931","NT5DC2"])

scatter(c3c4,"mEGFP-AP21967-2d","CDK2-AP21967-2d","SLAMseq.new.AP21967",highlight_genes=["CDK2","HBZ","ESPN","HBG2","HBG1","ENSG00000284931","C11orf21","TUBB1"])


# ----------------------------
df = pd.read_csv("SLAMseq.total.csv",header=0,index_col=0)
df = df.join(pd.read_csv("Gene.info.csv",index_col=0).iloc[:,0]).set_index("Symbol")
cpm = df.apply(lambda x:x*1000000/x.sum(),axis="index")

c1c2 = cpm[["C2.1","C2.2","C1.1","C1.2"]]
c3c4 = cpm[["C4.1","C4.2","C3.1","C3.2"]]

scatter(c1c2,"mEGFP-Rapa-2d","mEGFP-Rapa-2d","SLAMseq.total.Rapa",highlight_genes=["CDK2","HBZ","ESPN","HBG2","HBG1","ENSG00000284931","NT5DC2"])

scatter(c3c4,"mEGFP-AP21967-2d","CDK2-AP21967-2d","SLAMseq.total.AP21967",highlight_genes=["CDK2","HBZ","ESPN","HBG2","HBG1","ENSG00000284931","C11orf21","TUBB1"])