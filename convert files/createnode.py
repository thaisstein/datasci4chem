import xml.etree.ElementTree as ET 
import pandas as pd 
import os
import networkx as nx
import matplotlib.pyplot as plt
import csv
from matplotlib import pyplot as plt

def checkSmiles(smiles):
    if(smiles == " "):
        return 0
    return 1

def rowsLog(rows_log, regular, component, id, problematic_reactions):
    if(regular == 0):
        rows_log.append({"regular": regular,
        "reason": "nosmiles",
        "component" : component,
        "id_reac" : id
        })
        problematic_reactions.append(id)
#    else: 
#        rows_log.append({"regular": regular,
#        "reason": " ",
#        "component" : " ",
#        "id_reac" : " "
#        })


def rowsMolec(rows_molec, name, id):
    rows_molec.append({"id_molec": id, 
    "name": name.lower()
    }) 


def rowsReaction(rows_reac, smiles, id_reac, description, solvent, catalyst, problematic_reactions):
    if id_reac not in problematic_reactions:
        rows_reac.append({"id_reac": id_reac, 
                "smiles": smiles,
                "description": description, 
                "solvent": solvent, 
                "catalyst": catalyst
                }) 
            

def rowsEdge(source, target, role, rows_edge):
    rows_edge.append({"source": source,
    "target": target,
    "role" : role
    })
    
def rowsNewEdge(source, target, role, rows_edge):
    rows_edge.append({"source": source,
    "target": target,
    "id_reac" : role
    })

def csvFile(name, rows, cols):
    csvname = name  
    df = pd.DataFrame(rows, columns=cols)
    if(name.lower() == "moleculenode" or name.lower() =="reactionnode"):
        df = df.drop_duplicates()
    df.to_csv(csvname,  index = False, quoting=csv.QUOTE_ALL)

problematic_reactions = []

cols_reac = ["id_reac", "smiles", "description", "solvent", "catalyst", "temperature"] #TODO temperature
rows_reac = [] 

cols_molec = ["id_molec", "name"] 
rows_molec = []

cols_edge_01 = ["source", "target", "role"] 
rows_edge_01 = []

cols_edge_02 = ["source", "target", "role"] 
rows_edge_02 = []

cols_log = ["regular", "reason", "component", "id_reac"]
rows_log = []

cols_newedges = ["source", "target", "id"]
rows_newedges = []

l = ['pftaps19760106_wk01.xml' ,'pftaps19760113_wk02.xml' ,'pftaps19760120_wk03.xml', 'pftaps19760127_wk04.xml',
'pftaps19760203_wk05.xml', 'pftaps19760210_wk06.xml', 'pftaps19760217_wk07.xml', 'pftaps19760224_wk08.xml',
'pftaps19760302_wk09.xml', 'pftaps19760309_wk10.xml', 'pftaps19760316_wk11.xml', 'pftaps19760323_wk12.xml',
'pftaps19760330_wk13.xml', 'pftaps19760406_wk14.xml', 'pftaps19760413_wk15.xml', 'pftaps19760420_wk16.xml',
'pftaps19760427_wk17.xml', 'pftaps19760504_wk18.xml', 'pftaps19760511_wk19.xml', 'pftaps19760518_wk20.xml',
'pftaps19760525_wk21.xml', 'pftaps19760601_wk22.xml', 'pftaps19760608_wk23.xml', 'pftaps19760615_wk24.xml', 
'pftaps19760622_wk25.xml', 'pftaps19760629_wk26.xml', 'pftaps19760706_wk27.xml', 'pftaps19760713_wk28.xml',
'pftaps19760720_wk29.xml', 'pftaps19760727_wk30.xml', 'pftaps19760803_wk31.xml', 'pftaps19760810_wk32.xml',
'pftaps19760817_wk33.xml', 'pftaps19760824_wk34.xml', 'pftaps19760831_wk35.xml', 'pftaps19760907_wk36.xml',
'pftaps19760914_wk37.xml', 'pftaps19760921_wk38.xml', 'pftaps19760928_wk39.xml', 'pftaps19761005_wk40.xml',
'pftaps19761012_wk41.xml', 'pftaps19761019_wk42.xml', 'pftaps19761026_wk43.xml', 'pftaps19761102_wk44.xml',
'pftaps19761109_wk45.xml', 'pftaps19761116_wk46.xml', 'pftaps19761123_wk47.xml', 'pftaps19761130_wk48.xml',
'pftaps19761207_wk49.xml', 'pftaps19761214_wk50.xml', 'pftaps19761221_wk51.xml', 'pftaps19761228_wk52.xml']

for k in l:
    xmlfile = k # variable that keeps the file i'm using, will search for it in the same folder

    tree = ET.parse(xmlfile)
    root = tree.getroot() #reactionlist
    solvent = catalyst = " "


    id_reac = smiles = description = solvent = catalyst = temperature = id_reactants = id_products = name = source = target = role = " "

    count = {}  # dictionary for each patent
    for reaction in root: # each one of the reactions i the file
        productlist = []
        reactantlist = []
        for info in reaction:
            if (info.tag == "{http://bitbucket.org/dan2097}source"):
                for descript in info:
                    if(descript.tag == "{http://bitbucket.org/dan2097}documentId"):
                        id_reac = descript.text
                        if(id_reac in count): # adds id in dictionary
                            count[id_reac] += 1
                        else: 
                            count[id_reac] = 1 # id count starts with 1
                        id_reac = id_reac + str(count[id_reac])
                    if(descript.tag == "{http://bitbucket.org/dan2097}paragraphText"):
                        description = descript.text
            if (info.tag == "{http://bitbucket.org/dan2097}reactionSmiles"):
                smiles = info.text
            if (info.tag == "{http://www.xml-cml.org/schema}spectatorList"):
                for spectator in info:
                    if (spectator.attrib["role"] == "solvent"):
                        for i in spectator:
                            if(i.tag == "{http://www.xml-cml.org/schema}identifier" and i.attrib["dictRef"] == "cml:smiles"):
                                solvent = i.attrib["value"]

                    if (spectator.attrib["role"] == "catalyst"):
                        for j in spectator:
                            if(j.tag == "{http://www.xml-cml.org/schema}identifier" and j.attrib["dictRef"] == "cml:smiles"):
                                catalyst = j.text

            if (info.tag == "{http://www.xml-cml.org/schema}productList"): # products
                for product in info:
                    id_products = " "
                    smilesprod = 0
                    errorprod = 0
                    if (product.tag == "{http://www.xml-cml.org/schema}product"):

                        for prodinfo in product: # each product
                            if(prodinfo.tag == "{http://www.xml-cml.org/schema}molecule"):
                                for noun in prodinfo:
                                    if (noun.tag == "{http://bitbucket.org/dan2097}nameResolved"): 
                                        name = noun.text.lower()
                                    elif (noun.tag == "{http://www.xml-cml.org/schema}name"):
                                        name = noun.text.lower()
                            if (prodinfo.tag == "{http://www.xml-cml.org/schema}identifier" and prodinfo.attrib["dictRef"] == "cml:smiles"):
                                id_products = (prodinfo.attrib["value"])   
                                productlist.append(id_products)

                        smilesprod = checkSmiles(id_products)
                        rowsLog(rows_log, smilesprod, name, id_reac, problematic_reactions)
                        if(smilesprod!=0):
                            rowsMolec(rows_molec, name, id_products)
        
            if (info.tag == "{http://www.xml-cml.org/schema}reactantList"): # reactants
                for reactant in info:
                    id_reactants = " "
                    smilesreac = errorreac = 0
                    if (reactant.tag == "{http://www.xml-cml.org/schema}reactant"):
                        for reacinfo in reactant:
                            if(reacinfo.tag == "{http://www.xml-cml.org/schema}molecule"):
                                for noun in reacinfo:
                                    if (noun.tag == "{http://bitbucket.org/dan2097}nameResolved"): 
                                        name = noun.text.lower()
                                    elif (noun.tag == "{http://www.xml-cml.org/schema}name"):
                                        name = noun.text.lower()
                            if(reacinfo.tag == "{http://www.xml-cml.org/schema}identifier" and reacinfo.attrib["dictRef"]=="cml:smiles"):
                                id_reactants = (reacinfo.attrib["value"]) 
                                reactantlist.append(id_reactants)

                        smilesreac = checkSmiles(id_reactants)
                        rowsLog(rows_log, smilesreac, name, id_reac, problematic_reactions)
                        if (smilesreac != 0):
                            rowsMolec(rows_molec, name, id_reactants)    

        rowsReaction(rows_reac, smiles, id_reac, description, solvent, catalyst, problematic_reactions)    
       # print("product list", productlist)
        #print("reactant list", reactantlist)
       # print("prob", problematic_reactions)
        if(id_reac not in problematic_reactions):
            for j in range (len(productlist)):
                rowsEdge(id_reac, productlist[j], "product", rows_edge_01)
    
            for i in range (len(reactantlist)):
                rowsEdge(reactantlist[i], id_reac, "reactant", rows_edge_02)  
csvFile("moleculenode", rows_molec, cols_molec)
csvFile("reactionnode", rows_reac, cols_reac)
csvFile("edges_reactant_reaction", rows_edge_01, cols_edge_01)
csvFile("edges_reaction_product", rows_edge_02, cols_edge_02)
csvFile("log", rows_log, cols_log)


# erasing reactions on log


with open("log", 'r') as fh:
    total = 0
    irregular = 0
    id_list = []

    reader = csv.DictReader(fh, delimiter=',')

    for row in reader:
        total += 1
        reg = row['regular']
        id = row['id_reac']
        #print(reg)
        if(reg != '1'):
            irregular += 1
            id_list.append(id)

percentual = (irregular/total) * 100
format_percentual = "{:.2f}".format(percentual)

print(irregular, "of", total, "molecules are irregular")
print(format_percentual, "% of the molecules are invalid")

