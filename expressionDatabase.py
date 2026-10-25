import pandas as pd
import sqlite3
#Expression Atlas files from ebi.ac.uk, (gene, file path)
files = [("Lrrk", "expressionRawData/lrrk_expression.tsv"), ("Vps35", "expressionRawData/Vps35_expression.tsv"),("park", "expressionRawData/park_expression.tsv"), ("pink1", "expressionRawData/pink1_expression.tsv"), ("DJ1beta", "expressionRawData/DJ1beta_expression.tsv"), ("DJ1alpha", "expressionRawData/DJ1alpha_expression.tsv"), ("anne", "expressionRawData/anne_expression.tsv"),("iPLAVIA", "expressionRawData/iPLA2VIA_expression.tsv"), ("ntc", "expressionRawData/ntc_expression.tsv"), ("aux", "expressionRawData/aux_expression.tsv"), ("Synj", "expressionRawData/Synj_expression.tsv"), ("Vps13", "expressionRawData/Vps13_expression.tsv"), ("taf1", "expressionRawData/taf1_expression.tsv"),("Rab39","expressionRawData/Rab39_expression.tsv"),("RabX4","expressionRawData/RabX4_expression.tsv"), ("Gba1a", "expressionRawData/Gba1a_expression.tsv"), ("Gba1b", "expressionRawData/Gba1b_expression.tsv")]

table =[]

for gene, path in files:
    grid = pd.read_csv(path, sep="\t", comment ="#")
    grid = grid.rename(columns={grid.columns[0]: "experiment"})
    expression_long = grid.melt(id_vars = "experiment", var_name="cell_type", value_name="expression")
    expression_long["expression"] = pd.to_numeric(expression_long["expression"], errors = "coerce")
    experssion_long = expression_long.dropna(subset=["expression"])
    expression_long["gene"]=gene
    table.append(expression_long)

allGene = pd.concat(table,ignore_index =True)
group = allGene.groupby("gene")["expression"]
allGene["relative"]=(allGene["expression"] - group.transform("min")) / (group.transform("max") - group.transform("min"))
print(allGene.head())

conn = sqlite3.connect("pd_expression.db")
allGene.to_sql("expression", conn,if_exists = "replace", index=False)
print(pd.read_sql_query(
    "SELECT gene, cell_type, expression FROM expression "
    "WHERE cell_type LIKE '%PAM%' ORDER BY expression DESC", conn))

allGene.to_csv("expression_table.csv", index =False)
conn.close()


