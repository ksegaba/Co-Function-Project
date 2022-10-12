#!/usr/bin/env python3
"""
Description: Logistic Regression Model

Arguments:

Outputs:

"""
__author__ = "Kenia Segura Abá"

import sys, argparse
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix


if __name__=="__main__":
    # Argument Parser
    parser = argparse.ArgumentParser(description="")
    # Required input
    req_group = parser.add_argument_group(title="REQUIRED INPUT")
    req_group.add_argument("-x", help="Input feature matrix", required=True)
    req_group.add_argument("-test,", help="Test instances file", required=True)
    # Optional input
    opt_group = parser.add_argument_group(title="OPTIONAL INPUT")
    opt_group.add_argument("-y", help="Label file")
    opt_group.add_argument("-y_name", help="If not in separate file, give column name in file x")

