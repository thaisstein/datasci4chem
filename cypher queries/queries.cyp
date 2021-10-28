// loads molecules from moleculenode.csv
LOAD CSV WITH HEADERS FROM 'https://raw.githubusercontent.com/thaisstein/datasci4chemistry/main/moleculenode.csv' AS line
CREATE (:Molecule { id_molec: line.id_molec, name: line.name})

//creates index
CREATE INDEX FOR (n:Molecule) ON (n.code)


//loads reactions from reactionnode.csv

LOAD CSV WITH HEADERS FROM ‘https://raw.githubusercontent.com/thaisstein/datasci4chemistry/main/reactionnode.csv' AS line
CREATE (:Reaction { id_reac: line.id_reac, smiles: line.smiles, solvent: line.solvent, catalyst: line.catalyst})

//creates index
CREATE INDEX FOR (n:Reaction) ON (n.code)


//loads relation reaction to product
LOAD CSV WITH HEADERS FROM 'https://raw.githubusercontent.com/thaisstein/datasci4chemistry/main/edges_reaction_product.csv' AS line

MATCH (m:Molecule {id_molec: line.source})

MATCH (r:Reaction {id_reac: line.target})

CREATE (m)-[:Reacts {role: line.role}]->(r)


//loads relation reactant to reaction
LOAD CSV WITH HEADERS FROM 'https://raw.githubusercontent.com/thaisstein/datasci4chemistry/main/edges_reactant_reaction.csv' AS line

MATCH (r:Reaction {id_reac: line.source})
MATCH (m:Molecule {id_molec: line.target})


CREATE (r)-[:Generates {role: line.role}]->(m)



MATCH (d)-[:Reacts]->(p) RETURN d, p LIMIT 50


MATCH (m:Molecule)
WITH m.smiles AS smiles, COLLECT(m) AS branches
WHERE SIZE(branches) > 1
FOREACH (n IN branches | DETACH DELETE n);

//PROJECTION

MATCH(m1:Molecule)-[a]->(r1:Reaction)-[b]->(m2:Molecule)
WHERE m1.id_molec > m2.id_molec
WITH r1.id_reac AS reac
MERGE (m1)-[r:Relates{reaction: reac}]->(m2)

MATCH (m1:Molecule)<-[r:Relates]->(m2:Molecule)
RETURN m1, m2
LIMIT 50


