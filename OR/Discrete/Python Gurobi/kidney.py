import time
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from gurobipy import Model, GRB, quicksum

# ---------- Dataset ----------
pairs = [
    {"PairID": 1, "P_BT": "A", "P_PRA": 32, "P_HLA": "A1-B8-DR3", "P_CMV": "No", "P_Wt": 60, "P_Age": 45, "D_BT": "B", "D_HLA": "A2-B7-DR1", "D_CMV": "Yes", "D_Wt": 65, "D_Age": 42},
    {"PairID": 2, "P_BT": "B", "P_PRA": 77, "P_HLA": "A2-B15-DR4", "P_CMV": "Yes", "P_Wt": 55, "P_Age": 50, "D_BT": "A", "D_HLA": "A1-B8-DR3", "D_CMV": "No", "D_Wt": 65, "D_Age": 35},
    {"PairID": 3, "P_BT": "O", "P_PRA": 66, "P_HLA": "A3-B44-DR7", "P_CMV": "Yes", "P_Wt": 70, "P_Age": 33, "D_BT": "AB", "D_HLA": "A3-B44-DR1", "D_CMV": "Yes", "D_Wt": 80, "D_Age": 40},
    {"PairID": 4, "P_BT": "AB", "P_PRA": 28, "P_HLA": "A2-B18-DR2", "P_CMV": "No", "P_Wt": 65, "P_Age": 41, "D_BT": "O", "D_HLA": "A2-B18-DR1", "D_CMV": "No", "D_Wt": 60, "D_Age": 36},
    {"PairID": 5, "P_BT": "A", "P_PRA": 91, "P_HLA": "A1-B8-DR3", "P_CMV": "Yes", "P_Wt": 50, "P_Age": 38, "D_BT": "O", "D_HLA": "A2-B7-DR4", "D_CMV": "Yes", "D_Wt": 75, "D_Age": 50},
    {"PairID": 6, "P_BT": "B", "P_PRA": 35, "P_HLA": "A2-B15-DR5", "P_CMV": "No", "P_Wt": 68, "P_Age": 60, "D_BT": "AB", "D_HLA": "A1-B8-DR3", "D_CMV": "No", "D_Wt": 70, "D_Age": 34},
    {"PairID": 7, "P_BT": "O", "P_PRA": 85, "P_HLA": "A3-B7-DR2", "P_CMV": "Yes", "P_Wt": 60, "P_Age": 29, "D_BT": "A", "D_HLA": "A3-B7-DR2", "D_CMV": "Yes", "D_Wt": 64, "D_Age": 45},
    {"PairID": 8, "P_BT": "AB", "P_PRA": 42, "P_HLA": "A1-B8-DR1", "P_CMV": "No", "P_Wt": 59, "P_Age": 44, "D_BT": "O", "D_HLA": "A2-B18-DR4", "D_CMV": "Yes", "D_Wt": 61, "D_Age": 47},
    {"PairID": 9, "P_BT": "A", "P_PRA": 50, "P_HLA": "A2-B15-DR4", "P_CMV": "Yes", "P_Wt": 63, "P_Age": 52, "D_BT": "B", "D_HLA": "A2-B8-DR1", "D_CMV": "Yes", "D_Wt": 62, "D_Age": 46},
    {"PairID": 10, "P_BT": "B", "P_PRA": 38, "P_HLA": "A3-B7-DR3", "P_CMV": "No", "P_Wt": 58, "P_Age": 36, "D_BT": "A", "D_HLA": "A3-B7-DR1", "D_CMV": "No", "D_Wt": 55, "D_Age": 39},
    {"PairID": 11, "P_BT": "O", "P_PRA": 79, "P_HLA": "A1-B18-DR4", "P_CMV": "Yes", "P_Wt": 69, "P_Age": 31, "D_BT": "AB", "D_HLA": "A3-B7-DR2", "D_CMV": "Yes", "D_Wt": 72, "D_Age": 43},
    {"PairID": 12, "P_BT": "AB", "P_PRA": 45, "P_HLA": "A1-B7-DR3", "P_CMV": "No", "P_Wt": 62, "P_Age": 40, "D_BT": "B", "D_HLA": "A2-B15-DR2", "D_CMV": "No", "D_Wt": 67, "D_Age": 41},
    {"PairID": 13, "P_BT": "A", "P_PRA": 82, "P_HLA": "A3-B8-DR4", "P_CMV": "Yes", "P_Wt": 64, "P_Age": 48, "D_BT": "O", "D_HLA": "A2-B7-DR3", "D_CMV": "Yes", "D_Wt": 66, "D_Age": 49},
    {"PairID": 14, "P_BT": "B", "P_PRA": 31, "P_HLA": "A2-B15-DR1", "P_CMV": "No", "P_Wt": 61, "P_Age": 39, "D_BT": "A", "D_HLA": "A2-B8-DR1", "D_CMV": "No", "D_Wt": 64, "D_Age": 42},
    {"PairID": 15, "P_BT": "O", "P_PRA": 69, "P_HLA": "A3-B7-DR2", "P_CMV": "No", "P_Wt": 68, "P_Age": 34, "D_BT": "AB", "D_HLA": "A1-B18-DR3", "D_CMV": "Yes", "D_Wt": 70, "D_Age": 45},
    {"PairID": 16, "P_BT": "AB", "P_PRA": 55, "P_HLA": "A2-B44-DR4", "P_CMV": "Yes", "P_Wt": 66, "P_Age": 47, "D_BT": "O", "D_HLA": "A2-B44-DR2", "D_CMV": "No", "D_Wt": 68, "D_Age": 38},
    {"PairID": 17, "P_BT": "A", "P_PRA": 59, "P_HLA": "A1-B15-DR1", "P_CMV": "No", "P_Wt": 57, "P_Age": 53, "D_BT": "B", "D_HLA": "A2-B8-DR4", "D_CMV": "No", "D_Wt": 59, "D_Age": 51},
    {"PairID": 18, "P_BT": "B", "P_PRA": 41, "P_HLA": "A3-B7-DR5", "P_CMV": "Yes", "P_Wt": 60, "P_Age": 35, "D_BT": "A", "D_HLA": "A1-B15-DR2", "D_CMV": "Yes", "D_Wt": 62, "D_Age": 37},
    {"PairID": 19, "P_BT": "O", "P_PRA": 48, "P_HLA": "A2-B44-DR1", "P_CMV": "No", "P_Wt": 72, "P_Age": 30, "D_BT": "AB", "D_HLA": "A2-B7-DR3", "D_CMV": "Yes", "D_Wt": 74, "D_Age": 32},
    {"PairID": 20, "P_BT": "AB", "P_PRA": 29, "P_HLA": "A1-B8-DR2", "P_CMV": "Yes", "P_Wt": 67, "P_Age": 50, "D_BT": "O", "D_HLA": "A2-B8-DR4", "D_CMV": "No", "D_Wt": 69, "D_Age": 55},
    {"PairID": 21, "P_BT": 'AB', "P_PRA": 94, "P_HLA": 'B7-B8-A3', "P_CMV": 'Yes', "P_Wt": 62, "P_Age": 37, "D_BT": 'O', "D_HLA": 'B15-DR4-B44', "D_CMV": 'Yes', "D_Wt": 60, "D_Age": 59},
    {"PairID": 22, "P_BT": 'AB', "P_PRA": 84, "P_HLA": 'DR3-DR7-B18', "P_CMV": 'No', "P_Wt": 66, "P_Age": 27, "D_BT": 'B', "D_HLA": 'DR5-A3-DR4', "D_CMV": 'Yes', "D_Wt": 51, "D_Age": 37},
    {"PairID": 23, "P_BT": 'AB', "P_PRA": 32, "P_HLA": 'DR7-DR1-DR3', "P_CMV": 'Yes', "P_Wt": 68, "P_Age": 27, "D_BT": 'O', "D_HLA": 'DR4-A3-B7', "D_CMV": 'No', "D_Wt": 49, "D_Age": 55},
    {"PairID": 24, "P_BT": 'AB', "P_PRA": 37, "P_HLA": 'A1-DR3-B15', "P_CMV": 'No', "P_Wt": 69, "P_Age": 58, "D_BT": 'A', "D_HLA": 'B8-DR1-B44', "D_CMV": 'Yes', "D_Wt": 56, "D_Age": 60},
    {"PairID": 25, "P_BT": 'B', "P_PRA": 9, "P_HLA": 'DR2-A1-DR5', "P_CMV": 'No', "P_Wt": 68, "P_Age": 48, "D_BT": 'A', "D_HLA": 'DR4-A2-B7', "D_CMV": 'Yes', "D_Wt": 70, "D_Age": 42},
    {"PairID": 26, "P_BT": 'A', "P_PRA": 3, "P_HLA": 'DR2-DR5-DR1', "P_CMV": 'Yes', "P_Wt": 63, "P_Age": 40, "D_BT": 'O', "D_HLA": 'A1-DR7-DR5', "D_CMV": 'No', "D_Wt": 55, "D_Age": 53},
    {"PairID": 27, "P_BT": 'O', "P_PRA": 93, "P_HLA": 'DR2-A3-A2', "P_CMV": 'No', "P_Wt": 77, "P_Age": 33, "D_BT": 'O', "D_HLA": 'B18-DR5-B8', "D_CMV": 'No', "D_Wt": 55, "D_Age": 59},
    {"PairID": 28, "P_BT": 'A', "P_PRA": 35, "P_HLA": 'A3-DR4-B18', "P_CMV": 'Yes', "P_Wt": 63, "P_Age": 36, "D_BT": 'AB', "D_HLA": 'DR1-DR2-A1', "D_CMV": 'Yes', "D_Wt": 53, "D_Age": 41},
    {"PairID": 29, "P_BT": 'A', "P_PRA": 2, "P_HLA": 'DR4-DR2-DR3', "P_CMV": 'No', "P_Wt": 70, "P_Age": 28, "D_BT": 'A', "D_HLA": 'B7-DR7-B18', "D_CMV": 'Yes', "D_Wt": 47, "D_Age": 25},
    {"PairID": 30, "P_BT": 'AB', "P_PRA": 56, "P_HLA": 'DR2-DR5-B8', "P_CMV": 'Yes', "P_Wt": 59, "P_Age": 52, "D_BT": 'O', "D_HLA": 'B8-B44-DR7', "D_CMV": 'Yes', "D_Wt": 50, "D_Age": 51},
    {"PairID": 31, "P_BT": 'A', "P_PRA": 35, "P_HLA": 'A1-B18-B8', "P_CMV": 'No', "P_Wt": 62, "P_Age": 48, "D_BT": 'AB', "D_HLA": 'A2-DR3-DR1', "D_CMV": 'No', "D_Wt": 74, "D_Age": 27},
    {"PairID": 32, "P_BT": 'O', "P_PRA": 62, "P_HLA": 'DR7-B44-DR4', "P_CMV": 'No', "P_Wt": 70, "P_Age": 49, "D_BT": 'O', "D_HLA": 'DR3-DR4-A2', "D_CMV": 'Yes', "D_Wt": 68, "D_Age": 49},
    {"PairID": 33, "P_BT": 'A', "P_PRA": 10, "P_HLA": 'DR7-A2-B18', "P_CMV": 'No', "P_Wt": 64, "P_Age": 34, "D_BT": 'AB', "D_HLA": 'B15-B44-DR1', "D_CMV": 'No', "D_Wt": 61, "D_Age": 56},
    {"PairID": 34, "P_BT": 'A', "P_PRA": 100, "P_HLA": 'DR2-A3-DR4', "P_CMV": 'Yes', "P_Wt": 66, "P_Age": 53, "D_BT": 'B', "D_HLA": 'DR1-B44-DR2', "D_CMV": 'Yes', "D_Wt": 47, "D_Age": 31},
    {"PairID": 35, "P_BT": 'AB', "P_PRA": 38, "P_HLA": 'A1-DR5-B44', "P_CMV": 'Yes', "P_Wt": 57, "P_Age": 32, "D_BT": 'AB', "D_HLA": 'DR2-B44-A3', "D_CMV": 'Yes', "D_Wt": 51, "D_Age": 59},
    {"PairID": 36, "P_BT": 'B', "P_PRA": 80, "P_HLA": 'B15-A1-DR2', "P_CMV": 'Yes', "P_Wt": 80, "P_Age": 37, "D_BT": 'A', "D_HLA": 'DR5-A3-B7', "D_CMV": 'No', "D_Wt": 56, "D_Age": 29},
    {"PairID": 37, "P_BT": 'B', "P_PRA": 2, "P_HLA": 'B7-A3-A2', "P_CMV": 'No', "P_Wt": 62, "P_Age": 55, "D_BT": 'O', "D_HLA": 'DR3-A1-B15', "D_CMV": 'No', "D_Wt": 54, "D_Age": 44},
    {"PairID": 38, "P_BT": 'A', "P_PRA": 26, "P_HLA": 'B18-A3-B15', "P_CMV": 'Yes', "P_Wt": 76, "P_Age": 43, "D_BT": 'B', "D_HLA": 'B44-A2-DR4', "D_CMV": 'No', "D_Wt": 60, "D_Age": 53},
    {"PairID": 39, "P_BT": 'B', "P_PRA": 86, "P_HLA": 'DR4-DR3-A2', "P_CMV": 'Yes', "P_Wt": 64, "P_Age": 54, "D_BT": 'O', "D_HLA": 'A2-DR4-DR2', "D_CMV": 'Yes', "D_Wt": 56, "D_Age": 54},
    {"PairID": 40, "P_BT": 'B', "P_PRA": 31, "P_HLA": 'B15-B18-DR1', "P_CMV": 'No', "P_Wt": 59, "P_Age": 29, "D_BT": 'O', "D_HLA": 'A1-DR4-B8', "D_CMV": 'No', "D_Wt": 79, "D_Age": 45},
    {"PairID": 41, "P_BT": 'B', "P_PRA": 48, "P_HLA": 'B7-A3-A2', "P_CMV": 'Yes', "P_Wt": 79, "P_Age": 52, "D_BT": 'A', "D_HLA": 'DR4-DR2-A1', "D_CMV": 'Yes', "D_Wt": 51, "D_Age": 40},
    {"PairID": 42, "P_BT": 'B', "P_PRA": 10, "P_HLA": 'B44-DR4-A1', "P_CMV": 'Yes', "P_Wt": 63, "P_Age": 31, "D_BT": 'A', "D_HLA": 'B7-DR3-B15', "D_CMV": 'Yes', "D_Wt": 48, "D_Age": 41},
    {"PairID": 43, "P_BT": 'O', "P_PRA": 39, "P_HLA": 'B18-B44-DR1', "P_CMV": 'No', "P_Wt": 57, "P_Age": 51, "D_BT": 'A', "D_HLA": 'B15-B18-DR2', "D_CMV": 'No', "D_Wt": 80, "D_Age": 37},
    {"PairID": 44, "P_BT": 'A', "P_PRA": 51, "P_HLA": 'DR2-A1-B15', "P_CMV": 'Yes', "P_Wt": 64, "P_Age": 43, "D_BT": 'O', "D_HLA": 'A1-B8-DR4', "D_CMV": 'Yes', "D_Wt": 51, "D_Age": 31},
    {"PairID": 45, "P_BT": 'O', "P_PRA": 20, "P_HLA": 'DR5-B15-A3', "P_CMV": 'Yes', "P_Wt": 66, "P_Age": 47, "D_BT": 'O', "D_HLA": 'B8-DR2-B44', "D_CMV": 'Yes', "D_Wt": 46, "D_Age": 35},
    {"PairID": 46, "P_BT": 'O', "P_PRA": 18, "P_HLA": 'A1-DR1-B15', "P_CMV": 'No', "P_Wt": 51, "P_Age": 39, "D_BT": 'B', "D_HLA": 'A3-A2-B44', "D_CMV": 'Yes', "D_Wt": 59, "D_Age": 56},
    {"PairID": 47, "P_BT": 'A', "P_PRA": 90, "P_HLA": 'A3-A2-B7', "P_CMV": 'Yes', "P_Wt": 73, "P_Age": 50, "D_BT": 'B', "D_HLA": 'DR1-DR7-A2', "D_CMV": 'Yes', "D_Wt": 61, "D_Age": 58},
    {"PairID": 48, "P_BT": 'B', "P_PRA": 17, "P_HLA": 'DR3-A1-B8', "P_CMV": 'No', "P_Wt": 63, "P_Age": 54, "D_BT": 'O', "D_HLA": 'B15-A2-DR3', "D_CMV": 'Yes', "D_Wt": 62, "D_Age": 35},
    {"PairID": 49, "P_BT": 'A', "P_PRA": 17, "P_HLA": 'B7-A2-DR2', "P_CMV": 'Yes', "P_Wt": 53, "P_Age": 42, "D_BT": 'B', "D_HLA": 'DR3-DR5-DR1', "D_CMV": 'Yes', "D_Wt": 56, "D_Age": 29},
    {"PairID": 50, "P_BT": 'A', "P_PRA": 93, "P_HLA": 'DR1-B44-DR4', "P_CMV": 'No', "P_Wt": 54, "P_Age": 42, "D_BT": 'B', "D_HLA": 'A3-B8-B18', "D_CMV": 'Yes', "D_Wt": 62, "D_Age": 57},
    {"PairID": 51, "P_BT": 'B', "P_PRA": 34, "P_HLA": 'B15-B8-DR5', "P_CMV": 'Yes', "P_Wt": 60, "P_Age": 56, "D_BT": 'A', "D_HLA": 'B7-A2-DR7', "D_CMV": 'No', "D_Wt": 51, "D_Age": 53},
    {"PairID": 52, "P_BT": 'O', "P_PRA": 68, "P_HLA": 'DR5-B18-DR2', "P_CMV": 'No', "P_Wt": 76, "P_Age": 38, "D_BT": 'O', "D_HLA": 'B18-A2-B44', "D_CMV": 'Yes', "D_Wt": 65, "D_Age": 44},
    {"PairID": 53, "P_BT": 'B', "P_PRA": 8, "P_HLA": 'DR7-DR4-B15', "P_CMV": 'No', "P_Wt": 50, "P_Age": 38, "D_BT": 'AB', "D_HLA": 'DR2-A1-DR7', "D_CMV": 'Yes', "D_Wt": 52, "D_Age": 28},
    {"PairID": 54, "P_BT": 'O', "P_PRA": 87, "P_HLA": 'B18-B7-B8', "P_CMV": 'No', "P_Wt": 55, "P_Age": 29, "D_BT": 'O', "D_HLA": 'DR1-B7-A3', "D_CMV": 'Yes', "D_Wt": 78, "D_Age": 27},
    {"PairID": 55, "P_BT": 'AB', "P_PRA": 3, "P_HLA": 'DR2-DR1-DR7', "P_CMV": 'Yes', "P_Wt": 56, "P_Age": 40, "D_BT": 'AB', "D_HLA": 'B8-A1-DR3', "D_CMV": 'No', "D_Wt": 63, "D_Age": 56},
    {"PairID": 56, "P_BT": 'B', "P_PRA": 32, "P_HLA": 'B44-DR3-DR7', "P_CMV": 'Yes', "P_Wt": 69, "P_Age": 59, "D_BT": 'O', "D_HLA": 'DR4-DR2-B7', "D_CMV": 'No', "D_Wt": 51, "D_Age": 27},
    {"PairID": 57, "P_BT": 'O', "P_PRA": 92, "P_HLA": 'B8-DR3-B44', "P_CMV": 'No', "P_Wt": 72, "P_Age": 48, "D_BT": 'AB', "D_HLA": 'B15-DR7-B44', "D_CMV": 'No', "D_Wt": 83, "D_Age": 29},
    {"PairID": 58, "P_BT": 'A', "P_PRA": 41, "P_HLA": 'B15-DR5-DR1', "P_CMV": 'Yes', "P_Wt": 73, "P_Age": 49, "D_BT": 'AB', "D_HLA": 'B7-B44-A2', "D_CMV": 'No', "D_Wt": 78, "D_Age": 30},
    {"PairID": 59, "P_BT": 'B', "P_PRA": 49, "P_HLA": 'B18-A2-DR2', "P_CMV": 'Yes', "P_Wt": 56, "P_Age": 39, "D_BT": 'AB', "D_HLA": 'DR1-DR4-A2', "D_CMV": 'Yes', "D_Wt": 72, "D_Age": 25},
    {"PairID": 60, "P_BT": 'AB', "P_PRA": 81, "P_HLA": 'B44-DR7-A2', "P_CMV": 'No', "P_Wt": 70, "P_Age": 28, "D_BT": 'O', "D_HLA": 'B44-DR2-DR7', "D_CMV": 'No', "D_Wt": 85, "D_Age": 39},
    {"PairID": 61, "P_BT": 'A', "P_PRA": 50, "P_HLA": 'DR3-B8-B44', "P_CMV": 'Yes', "P_Wt": 80, "P_Age": 44, "D_BT": 'AB', "D_HLA": 'B18-A1-DR2', "D_CMV": 'Yes', "D_Wt": 79, "D_Age": 27},
    {"PairID": 62, "P_BT": 'O', "P_PRA": 42, "P_HLA": 'A3-B18-B7', "P_CMV": 'No', "P_Wt": 51, "P_Age": 57, "D_BT": 'B', "D_HLA": 'DR4-B18-B8', "D_CMV": 'No', "D_Wt": 54, "D_Age": 45},
    {"PairID": 63, "P_BT": 'A', "P_PRA": 57, "P_HLA": 'DR7-B44-A3', "P_CMV": 'Yes', "P_Wt": 68, "P_Age": 27, "D_BT": 'B', "D_HLA": 'DR1-A1-DR3', "D_CMV": 'Yes', "D_Wt": 57, "D_Age": 34},
    {"PairID": 64, "P_BT": 'O', "P_PRA": 12, "P_HLA": 'DR3-B18-A2', "P_CMV": 'Yes', "P_Wt": 58, "P_Age": 59, "D_BT": 'B', "D_HLA": 'DR7-A2-DR4', "D_CMV": 'Yes', "D_Wt": 47, "D_Age": 25},
    {"PairID": 65, "P_BT": 'O', "P_PRA": 40, "P_HLA": 'DR2-B7-DR5', "P_CMV": 'No', "P_Wt": 64, "P_Age": 60, "D_BT": 'A', "D_HLA": 'B7-DR3-DR7', "D_CMV": 'No', "D_Wt": 53, "D_Age": 31},
    {"PairID": 66, "P_BT": 'AB', "P_PRA": 35, "P_HLA": 'DR7-A2-B7', "P_CMV": 'No', "P_Wt": 53, "P_Age": 30, "D_BT": 'O', "D_HLA": 'DR1-B44-DR5', "D_CMV": 'No', "D_Wt": 71, "D_Age": 31},
    {"PairID": 67, "P_BT": 'B', "P_PRA": 78, "P_HLA": 'B44-B15-DR7', "P_CMV": 'Yes', "P_Wt": 67, "P_Age": 45, "D_BT": 'B', "D_HLA": 'A2-DR1-B7', "D_CMV": 'No', "D_Wt": 84, "D_Age": 44},
    {"PairID": 68, "P_BT": 'B', "P_PRA": 37, "P_HLA": 'B44-B7-DR4', "P_CMV": 'Yes', "P_Wt": 50, "P_Age": 33, "D_BT": 'A', "D_HLA": 'B18-DR2-A1', "D_CMV": 'Yes', "D_Wt": 64, "D_Age": 57},
    {"PairID": 69, "P_BT": 'AB', "P_PRA": 6, "P_HLA": 'DR4-A1-B7', "P_CMV": 'Yes', "P_Wt": 72, "P_Age": 55, "D_BT": 'B', "D_HLA": 'DR2-B18-DR7', "D_CMV": 'No', "D_Wt": 50, "D_Age": 55},
    {"PairID": 70, "P_BT": 'B', "P_PRA": 57, "P_HLA": 'A1-DR7-DR2', "P_CMV": 'No', "P_Wt": 79, "P_Age": 35, "D_BT": 'O', "D_HLA": 'DR4-B8-DR5', "D_CMV": 'No', "D_Wt": 85, "D_Age": 43},
    {"PairID": 71, "P_BT": 'AB', "P_PRA": 59, "P_HLA": 'DR7-B7-A2', "P_CMV": 'Yes', "P_Wt": 67, "P_Age": 60, "D_BT": 'B', "D_HLA": 'A3-DR2-DR5', "D_CMV": 'No', "D_Wt": 56, "D_Age": 54},
    {"PairID": 72, "P_BT": 'A', "P_PRA": 9, "P_HLA": 'A1-DR4-DR2', "P_CMV": 'No', "P_Wt": 62, "P_Age": 51, "D_BT": 'B', "D_HLA": 'DR5-B18-B8', "D_CMV": 'Yes', "D_Wt": 48, "D_Age": 40},
    {"PairID": 73, "P_BT": 'AB', "P_PRA": 58, "P_HLA": 'A1-DR7-A3', "P_CMV": 'No', "P_Wt": 64, "P_Age": 45, "D_BT": 'B', "D_HLA": 'DR5-B7-B18', "D_CMV": 'Yes', "D_Wt": 59, "D_Age": 30},
    {"PairID": 74, "P_BT": 'B', "P_PRA": 89, "P_HLA": 'A1-DR3-B44', "P_CMV": 'No', "P_Wt": 68, "P_Age": 26, "D_BT": 'B', "D_HLA": 'DR5-DR4-DR3', "D_CMV": 'No', "D_Wt": 54, "D_Age": 33},
    {"PairID": 75, "P_BT": 'O', "P_PRA": 14, "P_HLA": 'A3-B18-DR7', "P_CMV": 'Yes', "P_Wt": 55, "P_Age": 55, "D_BT": 'B', "D_HLA": 'DR1-DR3-DR2', "D_CMV": 'Yes', "D_Wt": 78, "D_Age": 56},
    {"PairID": 76, "P_BT": 'B', "P_PRA": 91, "P_HLA": 'B8-DR3-B15', "P_CMV": 'Yes', "P_Wt": 72, "P_Age": 45, "D_BT": 'A', "D_HLA": 'B44-B15-A3', "D_CMV": 'Yes', "D_Wt": 75, "D_Age": 60},
    {"PairID": 77, "P_BT": 'O', "P_PRA": 44, "P_HLA": 'A1-B15-B8', "P_CMV": 'No', "P_Wt": 56, "P_Age": 37, "D_BT": 'AB', "D_HLA": 'DR7-DR4-A1', "D_CMV": 'Yes', "D_Wt": 46, "D_Age": 46},
    {"PairID": 78, "P_BT": 'A', "P_PRA": 10, "P_HLA": 'B15-DR2-DR5', "P_CMV": 'No', "P_Wt": 63, "P_Age": 27, "D_BT": 'B', "D_HLA": 'B7-DR2-DR7', "D_CMV": 'Yes', "D_Wt": 67, "D_Age": 42},
    {"PairID": 79, "P_BT": 'A', "P_PRA": 92, "P_HLA": 'B7-B15-DR5', "P_CMV": 'Yes', "P_Wt": 58, "P_Age": 31, "D_BT": 'O', "D_HLA": 'B15-DR3-DR1', "D_CMV": 'Yes', "D_Wt": 51, "D_Age": 39},
    {"PairID": 80, "P_BT": 'A', "P_PRA": 1, "P_HLA": 'B18-DR1-DR7', "P_CMV": 'Yes', "P_Wt": 59, "P_Age": 44, "D_BT": 'AB', "D_HLA": 'A3-B7-A1', "D_CMV": 'No', "D_Wt": 48, "D_Age": 26},
    {"PairID": 81, "P_BT": 'B', "P_PRA": 94, "P_HLA": 'DR5-B7-B8', "P_CMV": 'Yes', "P_Wt": 59, "P_Age": 29, "D_BT": 'AB', "D_HLA": 'B44-A3-B18', "D_CMV": 'Yes', "D_Wt": 83, "D_Age": 39},
    {"PairID": 82, "P_BT": 'A', "P_PRA": 80, "P_HLA": 'A3-B8-A2', "P_CMV": 'Yes', "P_Wt": 65, "P_Age": 46, "D_BT": 'O', "D_HLA": 'DR2-B7-A1', "D_CMV": 'Yes', "D_Wt": 68, "D_Age": 40},
    {"PairID": 83, "P_BT": 'O', "P_PRA": 6, "P_HLA": 'B8-DR3-B18', "P_CMV": 'No', "P_Wt": 69, "P_Age": 45, "D_BT": 'B', "D_HLA": 'A3-B15-DR4', "D_CMV": 'Yes', "D_Wt": 66, "D_Age": 58},
    {"PairID": 84, "P_BT": 'B', "P_PRA": 25, "P_HLA": 'A1-B15-DR1', "P_CMV": 'Yes', "P_Wt": 52, "P_Age": 37, "D_BT": 'AB', "D_HLA": 'DR2-DR1-A3', "D_CMV": 'No', "D_Wt": 61, "D_Age": 57},
    {"PairID": 85, "P_BT": 'O', "P_PRA": 65, "P_HLA": 'B8-DR5-B44', "P_CMV": 'No', "P_Wt": 72, "P_Age": 32, "D_BT": 'A', "D_HLA": 'A3-A2-DR1', "D_CMV": 'No', "D_Wt": 73, "D_Age": 28},
    {"PairID": 86, "P_BT": 'O', "P_PRA": 29, "P_HLA": 'B44-DR3-A3', "P_CMV": 'No', "P_Wt": 61, "P_Age": 32, "D_BT": 'A', "D_HLA": 'DR7-B15-A1', "D_CMV": 'Yes', "D_Wt": 55, "D_Age": 50},
    {"PairID": 87, "P_BT": 'O', "P_PRA": 43, "P_HLA": 'DR5-DR3-B8', "P_CMV": 'Yes', "P_Wt": 80, "P_Age": 57, "D_BT": 'AB', "D_HLA": 'B7-DR7-B15', "D_CMV": 'No', "D_Wt": 75, "D_Age": 58},
    {"PairID": 88, "P_BT": 'B', "P_PRA": 91, "P_HLA": 'DR4-DR5-A3', "P_CMV": 'No', "P_Wt": 75, "P_Age": 25, "D_BT": 'A', "D_HLA": 'DR7-A2-DR3', "D_CMV": 'Yes', "D_Wt": 56, "D_Age": 29},
    {"PairID": 89, "P_BT": 'A', "P_PRA": 50, "P_HLA": 'DR3-DR5-B15', "P_CMV": 'No', "P_Wt": 79, "P_Age": 37, "D_BT": 'O', "D_HLA": 'DR4-B15-A3', "D_CMV": 'Yes', "D_Wt": 50, "D_Age": 40},
    {"PairID": 90, "P_BT": 'A', "P_PRA": 84, "P_HLA": 'A1-B18-A3', "P_CMV": 'Yes', "P_Wt": 73, "P_Age": 50, "D_BT": 'AB', "D_HLA": 'A2-A3-B15', "D_CMV": 'No', "D_Wt": 69, "D_Age": 27},
    {"PairID": 91, "P_BT": 'A', "P_PRA": 56, "P_HLA": 'DR5-B44-B8', "P_CMV": 'No', "P_Wt": 61, "P_Age": 42, "D_BT": 'B', "D_HLA": 'B18-B7-B8', "D_CMV": 'No', "D_Wt": 78, "D_Age": 55},
    {"PairID": 92, "P_BT": 'AB', "P_PRA": 32, "P_HLA": 'DR1-A1-DR2', "P_CMV": 'Yes', "P_Wt": 72, "P_Age": 33, "D_BT": 'O', "D_HLA": 'DR2-B15-B7', "D_CMV": 'No', "D_Wt": 45, "D_Age": 34},
    {"PairID": 93, "P_BT": 'AB', "P_PRA": 52, "P_HLA": 'A2-B7-B15', "P_CMV": 'No', "P_Wt": 60, "P_Age": 45, "D_BT": 'B', "D_HLA": 'DR3-A3-A1', "D_CMV": 'No', "D_Wt": 69, "D_Age": 53},
    {"PairID": 94, "P_BT": 'A', "P_PRA": 26, "P_HLA": 'A1-B18-DR5', "P_CMV": 'Yes', "P_Wt": 73, "P_Age": 26, "D_BT": 'B', "D_HLA": 'B8-DR7-A1', "D_CMV": 'No', "D_Wt": 53, "D_Age": 50},
    {"PairID": 95, "P_BT": 'O', "P_PRA": 37, "P_HLA": 'A2-DR4-DR5', "P_CMV": 'Yes', "P_Wt": 72, "P_Age": 43, "D_BT": 'A', "D_HLA": 'B7-B44-B18', "D_CMV": 'Yes', "D_Wt": 50, "D_Age": 35},
    {"PairID": 96, "P_BT": 'B', "P_PRA": 95, "P_HLA": 'DR5-DR4-DR1', "P_CMV": 'No', "P_Wt": 50, "P_Age": 39, "D_BT": 'A', "D_HLA": 'DR4-DR7-A2', "D_CMV": 'Yes', "D_Wt": 67, "D_Age": 54},
    {"PairID": 97, "P_BT": 'AB', "P_PRA": 0, "P_HLA": 'A3-B7-DR2', "P_CMV": 'Yes', "P_Wt": 64, "P_Age": 53, "D_BT": 'B', "D_HLA": 'B15-B44-A2', "D_CMV": 'Yes', "D_Wt": 45, "D_Age": 31},
    {"PairID": 98, "P_BT": 'B', "P_PRA": 5, "P_HLA": 'DR4-B7-A1', "P_CMV": 'Yes', "P_Wt": 55, "P_Age": 53, "D_BT": 'B', "D_HLA": 'B8-DR2-B15', "D_CMV": 'Yes', "D_Wt": 72, "D_Age": 29},
    {"PairID": 99, "P_BT": 'B', "P_PRA": 49, "P_HLA": 'B8-A1-DR7', "P_CMV": 'No', "P_Wt": 71, "P_Age": 55, "D_BT": 'B', "D_HLA": 'DR5-A1-B8', "D_CMV": 'No', "D_Wt": 54, "D_Age": 32},
    {"PairID": 100, "P_BT": 'A', "P_PRA": 80, "P_HLA": 'A1-A3-DR2', "P_CMV": 'No', "P_Wt": 66, "P_Age": 37, "D_BT": 'O', "D_HLA": 'B7-B18-A1', "D_CMV": 'Yes', "D_Wt": 51, "D_Age": 49}
]


df_pairs = pd.DataFrame(pairs)

donors = [
    {'DonorID': 103, 'D_BT': 'B', 'D_HLA': 'A1-B28-DR1', 'D_CMV': 'Yes', 'D_Wt': 64, 'D_Age': 34},
    {'DonorID': 109, 'D_BT': 'AB', 'D_HLA': 'A1-B12-DR2', 'D_CMV': 'No', 'D_Wt': 73, 'D_Age': 23},
    {'DonorID': 110, 'D_BT': 'O', 'D_HLA': 'A2-B18-DR4', 'D_CMV': 'No', 'D_Wt': 81, 'D_Age': 53}   
]

df_donors = pd.DataFrame(donors)

# ---------- Compatibility Rules ----------
bt_compat = {
    "O": ["O", "A", "B", "AB"],
    "A": ["A", "AB"],
    "B": ["B", "AB"],
    "AB": ["AB"]
}

def hla_mismatch(hla1, hla2):
    """
    Calculate the HLA mismatch score between two HLA strings.

    The HLA strings are expected to be in the format like "A1-B8-DR3",
    where the components are separated by dashes.

    Parameters:
    -----------
    hla1 : str
        HLA string for the first subject (e.g., donor).
    hla2 : str
        HLA string for the second subject (e.g., patient).

    Returns:
    --------
    int
        The number of mismatched HLA components, computed as the size
        of the symmetric difference between the sets of components.
    """
    return len(set(hla1.split("-")).symmetric_difference(set(hla2.split("-"))))

def is_compatible(donor, patient):
    """
    Determine if a donor is compatible with a patient based on
    blood type, PRA level, CMV status, weight, and HLA mismatch.

    Compatibility rules:
    - Patient's blood type must be compatible with donor's blood type.
    - Patient's PRA (panel reactive antibody) must be <= 80.
    - Donor with positive CMV (Cytomegalovirus) cannot donate to CMV-negative patient.
    - Donor's weight must be >= patient's weight.
    - HLA mismatch score must be 3 or less.

    Parameters:
    -----------
    donor : dict
        Dictionary containing donor attributes with keys:
        'D_BT', 'D_HLA', 'D_CMV', 'D_Wt', etc.
    patient : dict
        Dictionary containing patient attributes with keys:
        'P_BT', 'P_PRA', 'P_HLA', 'P_CMV', 'P_Wt', etc.

    Returns:
    --------
    bool
        True if donor is compatible with patient, False otherwise.
    """
    # TODO
    return True

# ---------- Build Graph ----------
def build_graph(df_pairs, df_donors=pd.DataFrame(), pruning= True):
    """
    Build a directed compatibility graph from patient-donor pairs and
    optionally altruistic donors.

    Nodes represent patient-donor pairs and altruistic donors.
    Directed edges represent compatibility from a donor to a patient.

    Parameters:
    -----------
    df_pairs : pandas.DataFrame
        DataFrame containing patient-donor pairs with at least columns:
        'PairID', donor and patient attributes.
    df_donors : pandas.DataFrame, optional
        DataFrame containing altruistic donors with donor attributes
        (default is an empty DataFrame).
    pruning : bool, optional
        Whether to prune the graph by removing incompatible nodes
        (default is True).

    Returns:
    --------
    networkx.DiGraph
        A directed graph with nodes as pairs/donors and edges indicating
        donor-to-patient compatibility.
    """
    graph = nx.DiGraph()
    # Add a node for each patient-donor pair
    for idx in df_pairs['PairID']:
        graph.add_node(idx)

    compatible_edges = []
    # Add edges between patient-donor pairs
    for donor in df_pairs.itertuples():
        for patient in df_pairs.itertuples():
            if donor.PairID == patient.PairID:
                continue
            if is_compatible(donor._asdict(), patient._asdict()):
                compatible_edges.append((donor.PairID, patient.PairID))
                graph.add_edge(donor.PairID, patient.PairID)

    if not df_donors.empty:
        # Add nodes for altruistic donors to the graph
        for idx in df_donors['DonorID']:
            graph.add_node(idx)

        # Add edges from altruistic donors to compatible patients
        for donor in df_donors.itertuples():
            for patient in df_pairs.itertuples():
                if is_compatible(donor._asdict(), patient._asdict()):
                    compatible_edges.append((donor.DonorID, patient.PairID))
                    graph.add_edge(donor.DonorID, patient.PairID)

    if pruning:
        prune_nodes(graph, df_donors)

    return graph

def prune_nodes(graph: nx.DiGraph, df_donors=pd.DataFrame()):
    """
    Recursively remove nodes from the graph that have either zero in-degree
    or zero out-degree, except:
      - Altruistic Donors that have outgoing edges are preserved.
      - Nodes reachable from altruistic donors are preserved.

    This pruning aims to remove isolated or dead-end nodes, improving
    the structure of the compatibility graph.

    Parameters:
    -----------
    graph : networkx.DiGraph
        The directed compatibility graph to prune.
    df_donors : pandas.DataFrame, optional
        DataFrame containing altruistic donors with donor attributes
        (default is an empty DataFrame).

    Returns:
    --------
    None
        The graph is modified in-place.
    """
    # TODO

def plot_graph(graph: nx.DiGraph, highlight_edges=None):
    """
    Visualize the donor-patient compatibility graph using matplotlib.

    Nodes are shown with labels, and edges represent compatibility.

    Parameters:
    -----------
    graph : networkx.DiGraph
        The compatibility graph to visualize.
    highlight_edges: list
        The list of edges to be highlighted.
    Returns:
    --------
    None
    """
    plt.figure(figsize=(10, 7))
    pos = nx.kamada_kawai_layout(graph)

    # Draw all nodes
    nx.draw_networkx_nodes(graph, pos, node_color='lightgreen', node_size=500)
    nx.draw_networkx_labels(graph, pos)

    # Separate edges into highlighted and normal
    highlight_edges = set(highlight_edges) if highlight_edges else set()
    all_edges = set(graph.edges)
    normal_edges = all_edges - highlight_edges

    # Draw normal edges in black
    nx.draw_networkx_edges(graph, pos, edgelist=list(normal_edges), edge_color='black', arrows=True, arrowsize=15)

    # Draw highlighted edges in red
    if highlight_edges:
        nx.draw_networkx_edges(graph, pos, edgelist=list(highlight_edges), edge_color='red', arrows=True, arrowsize=15, width=2)

    plt.title("Donor → Patient Compatibility Graph")
    plt.axis('off')
    plt.show()

# -------------------------------------
def recursive_algorithm(graph: nx.DiGraph, df_pairs, df_donors, max_cycle):
    """
    Solves a kidney exchange optimization problem using integer programming.

    The function models and solves a directed kidney exchange graph, selecting the
    optimal set of edges (i.e., transplants) that maximize the number of transplants
    while obeying flow constraints and avoiding long cycles.

    Parameters:
    -----------
    graph : nx.DiGraph
        A directed graph where nodes represent donor-patient pairs or altruistic donors,
        and edges represent feasible transplantations.

    df_pairs : pandas.DataFrame
        DataFrame containing paired donor-patient information.
        Must contain a 'PairID' column matching nodes in the graph.

    df_donors : pandas.DataFrame
        DataFrame containing altruistic donor information.
        Must contain a 'DonorID' column matching nodes in the graph.

    max_cycle : int
        The maximum allowed cycle length in the solution. Longer cycles
        are iteratively removed by adding constraints.

    Returns:
    --------
    selected_edges : list of tuple
        A list of directed edges (tuples) representing the selected transplant paths.
        If the model is infeasible or no optimal solution is found, returns `None`.

    Notes:
    ------
    - Uses Gurobi to solve a mixed-integer program.
    - Flow variables (`f_in`, `f_out`) are used to ensure logical consistency in edge selection.
    - Constraints:
        - Each node has balanced in/out flow.
        - Pair nodes may have one incoming and one outgoing edge.
        - Donor nodes may have at most one outgoing edge.
    - Iteratively detects and removes cycles of length ≥ `max_cycle`.
    """
    all_nodes = set(graph.nodes)
    
    # Filter df_pairs and df_donors to match the current graph
    df_pairs = df_pairs[df_pairs['PairID'].isin(all_nodes)].copy()
    df_donors = df_donors[df_donors['DonorID'].isin(all_nodes)].copy()

    pair_nodes = set(df_pairs['PairID'].values)
    donor_nodes = set(df_donors['DonorID'].values)
    

    model = Model("KidneyExchange")
    model.setParam("OutputFlag", 0)  # Optional: silence solver output

    # y_e binary variables for each edge
    y = {
        e: model.addVar(vtype=GRB.BINARY, name=f"y_{e}")
        for e in graph.edges
    }

    # Flow indicator variables per node
    f_in = {v: model.addVar(name=f"f_in_{v}") for v in all_nodes}
    f_out = {v: model.addVar(name=f"f_out_{v}") for v in all_nodes}

    model.update()

    # Link flow indicators to edge usage
    for v in all_nodes:
        in_edges = list(graph.in_edges(v))
        out_edges = list(graph.out_edges(v))

        model.addConstr(quicksum(y[e] for e in in_edges) == f_in[v])
        model.addConstr(quicksum(y[e] for e in out_edges) == f_out[v])

    # Pair node flow constraints
    for v in pair_nodes:
        model.addConstr(f_out[v] <= f_in[v])
        model.addConstr(f_in[v] <= 1)

    # Donor node flow constraints
    for v in donor_nodes:
        model.addConstr(f_out[v] <= 1)

    # Objective: maximize number of selected edges
    model.setObjective(quicksum(y.values()), GRB.MAXIMIZE)

    # Iteratively eliminate long cycles
    cycles_eliminated = 0
    while True:
        model.optimize()

        if model.status != GRB.OPTIMAL:
            print("No optimal solution found.")
            return None

        selected_edges = [e for e in y if y[e].X > 0.5]

        # Check for long cycles
        subgraph = nx.DiGraph()
        subgraph.add_edges_from(selected_edges)
        cycles = list(nx.simple_cycles(subgraph))
        long_cycles = [c for c in cycles if len(c) > max_cycle]

        if not long_cycles:
            print("Number of removed Cycles: {}".format(cycles_eliminated))
            return selected_edges  # Final solution

        # TODO: Remove all cycles in long_cycles
        model.update()

# ------------------------------------------------------------

if __name__ == "__main__":
    G = build_graph(df_pairs, pruning=False)
    print(G.number_of_edges())
    G = build_graph(df_pairs, pruning=True)
    print(G.number_of_edges())
    solution = recursive_algorithm(G, df_pairs, df_donors, 7)
    print("Number of kidneys transplanted: {}".format(len(solution)))
    solution = recursive_algorithm(G, df_pairs, df_donors, 3)
    print("Number of kidneys transplanted: {}".format(len(solution)))
    G = build_graph(df_pairs, df_donors, pruning=True)
    solution = recursive_algorithm(G, df_pairs, df_donors, 3)
    print("Number of kidneys transplanted: {}".format(len(solution)))
    plot_graph(G, solution)