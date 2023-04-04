"""
Pre-process SGD reaction data in the biopax format
"""

import os
from bs4 import BeautifulSoup

os.chdir("/mnt/home/seguraab/Shiu_Lab/Co-function/")

class Reaction:
    def __init__(self, head_tag):
        """
        Instantiate an instance of the Reaction class.
        Parameters:
            gene, str: gene ID of enzyme catalyzing the reaction
            protein, str: protein standard name
            direction, str: direction of reaction
            step_num, int: indicates the pathway step number when the reaction occurs
            pathway, str: ID of the pathway the reaction is a part of
        Returns: none
        """
        self.head_tag = head_tag
        self.reactants = [] # list of reactants in the reaction
        self.products = [] # list of products in the reaction
        self.rxn_ids = {}
        self.direction = None
        self.step_num = None # idk how to incorporate this information (maybe use the step ID and then at the end assign step num based on the Pathway leaf?)


class Pathway:
    """
    Class made up of a series of reactions
    """



def parse_BiochemicalPathwayStep(child):
    """
    Parse the BiochemicalPathwayStep child in the biopax pathway tree
    return None
    """
    try:
        for descendant in child.children:
            # get the ID of the reaction
            if descendant.name == "stepConversion":
                step_id = descendant["rdf:resource"]
                rxn_ids[step_id] = {
                    "rxn":{},
                    "neighbors":[] # or should it be a dict of step_ids
                }
            # get the direction of the reaction
            if descendant.name == "stepDirection": # direction of step
                print(f"Descendant name: {descendant.name}")
                rxn_ids[step_id]["direction"] = descendant.get_text()
            # get gene name of protein catalyzing the reaction
            elif descendant.name == "stepProcess":
                #rxn = Reaction() # instantiate a reaction
                # print(descendant.controller.Protein.prettify())
                rxn_id = descendant.Catalysis["rdf:ID"]
                gene = descendant.controller.Protein.xref.id.get_text()
                protein = descendant.controller.Protein.standardName.get_text()
                rxn_ids[step_id]["rxn"][rxn_id] = [gene, protein]
    except TypeError:
        print(f"This child is in the wrong format:\n{child.prettify()}")


def parse_pathway(child, rxn_ids):
    """
    Parse the pathway steps in the biopax pathway tree
    """
    try:
        for descendant in child.children:
            if child.name == "BiochemicalPathwayStep": # for each step in pathway
                parse_BiochemicalPathwayStep(child)
            # determine if the current step has a neighboring reaction
            elif descendant.name == "nextStep":
                print(descendant.BiochemicalPathwayStep.name)
                parse_BiochemicalPathwayStep(descendant)
                neighbor_id = descendant.BiochemicalPathwayStep["rdf:ID"]
                rxn_ids["neighbors"].append(neighbor_id)
    except TypeError:
        print(f"This child is in the wrong format:\n{child.prettify()}")


        

if __name__ == "__main__":
    # Read XML files (args or a loop?)
    with open("Data/SGD_BioPax/GLYCOLYSIS.xml", "r") as f:
        data = f.read()
        soup = BeautifulSoup(data, "xml")
    
    head_tag = soup.RDF # tree head
    rxn_ids = {} # reactions in the pathway
    
    # parse_pathway(child, rxn_ids)
    for child in head_tag.children:
        print(f"Child name: {child.name}")
        if child.name == "BiochemicalPathwayStep": # for each step in pathway
            parse_BiochemicalPathwayStep (child)
            # determine if the current step has a neighboring reaction
            for descendant in child.children:
                if descendant.name == "nextStep":
                    print(descendant.name)
                    # print(descendant.BiochemicalPathwayStep)
                    parse_BiochemicalPathwayStep(descendant)
                    neighbor_id = descendant.BiochemicalPathwayStep["rdf:ID"]
                    rxn_ids["neighbors"].append(neighbor_id)
