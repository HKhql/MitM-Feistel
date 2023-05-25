

from utils import *
from gurobipy import *
import os
import sys


class Vars_generator:
    @staticmethod
    def genVars_input_of_round(i, r, pos):
        return [f'IP_{i}_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_output_of_round(i, r, pos):
        return [f'OP_{i}_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_output_of_ShiftRow(i, r, pos):
        return [f'OSR_{i}_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_output_of_Xor(i, r, pos):
        return [f'OXor_{i}_r{r}_{pos}_{j}' for j in range(bs)]

    # Separate
    @staticmethod
    def genVars_In_isWhite(r, op):
        return [f'In_isWhite_r{r}_{op}_{j}' for j in range(bs)]

    # MixColumn
    @staticmethod
    def genVars_SupP_Blue_input_of_MixColumn(i, r, op):
        return [f'IMC_SupP_Blue_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Red_input_of_MixColumn(i, r, op):
        return [f'IMC_SupP_Red_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Blue_output_of_MixColumn(i, r, op):
        return [f'OMC_SupP_Blue_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Red_output_of_MixColumn(i, r, op):
        return [f'OMC_SupP_Red_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_MC_SupP_Blue_ColExistWhite(r, op):
        return [f'MC_SupP_Blue_ColExistWhite_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Red_ColExistWhite(r, op):
        return [f'MC_SupP_Red_ColExistWhite_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Blue_ColAllGray(r, op):
        return [f'MC_SupP_Blue_ColAllGray_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Red_ColAllGray(r, op):
        return [f'MC_SupP_Red_ColAllGray_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Blue_SumGray(r, op):
        return [f'G_SupP_Blue_SumGray_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Red_SumGray(r, op):
        return [f'G_SupP_Red_SumGray_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_ConsumedDeg_Blue(r, op):
        return [f'G_CD_MC_Blue_r{r}_{op}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_ConsumedDeg_Red(r, op):
        return [f'G_CD_MC_Red__r{r}_{op}_{j}' for j in range(ColN)]

    # Xor
    @staticmethod
    def genVars_Xor_ConsumedDeg_Blue(r, op):
        return [f'CD_Xor_Blue_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Xor_ConsumedDeg_Red(r, op):
        return [f'CD_Xor_Red__r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Blue_input_of_Xor(i, r, op):
        return [f'IXor_SupP_Blue_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Red_input_of_Xor(i, r, op):
        return [f'IXor_SupP_Red_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Blue_output_of_Xor(i, r, op):
        return [f'OXor_SupP_Blue_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Red_output_of_Xor(i, r, op):
        return [f'OXor_SupP_Red_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_isWhite(r, op):
        return [f'OXor_isWhite_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Blue_AND(i, r, op):
        return [f'OXor_SupP_Blue_AND_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Blue_OR(i, r, op):
        return [f'OXor_SupP_Blue_OR_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Red_AND(i, r, op):
        return [f'OXor_SupP_Red_AND_{i}_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Red_OR(i, r, op):
        return [f'OXor_SupP_Red_OR_{i}_r{r}_{op}_{j}' for j in range(bs)]

    # Initial Degree
    @staticmethod
    def genVars_degree__forward():
        return ['deg_f_' + str(j) for j in range(bs * b)]

    @staticmethod
    def genVars_degree_backward():
        return ['deg_b_' + str(j) for j in range(bs * b)]

    # Additional
    @staticmethod
    def genVasr_number_of_initial_degree():
        return 'Indicator_degree_less'

    # Match
    @staticmethod
    def genVar_Match_Counted(pos):
        return [f'G_Match_Counted_{pos}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_Match_Counted_direct():
        return [f'Match_Counted_direct_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_IMC_isWhite(pos):
        return [f'Match_IMC_isWhite_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_IP_isWhite(pos):
        return [f'Match_IP_isWhite_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_OP_isWhite(pos):
        return [f'Match_OP_isWhite_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_OXor_isWhite(pos):
        return [f'Match_OXor_isWhite_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_Exist(pos):
        return [f'Match_Exist_{pos}_{j}' for j in range(ColN)]


class Constraints_generator:
    def __init__(self, total_round, initial_round, matching_round):
        self.ini_r = initial_round
        self.mat_r = matching_round
        self.TR = total_round
        self.pos4 = ['LL', 'LR', 'RL', 'RR']
        self.pos3 = ['L', 'R', 'M']
        self.op = ['MCL', 'MCR', 'XorL', 'XorR', 'MCM']
        self.output_index = output_index(self.TR)

    def genConstraints_initial_degree(self):
        cons = []
        IP_1 = []
        IP_2 = []
        d1 = Vars_generator.genVars_degree__forward()
        d2 = Vars_generator.genVars_degree_backward()
        for pos in self.pos4:
            IP_1 = IP_1 + Vars_generator.genVars_input_of_round(1, self.ini_r, pos)
            IP_2 = IP_2 + Vars_generator.genVars_input_of_round(2, self.ini_r, pos)
        for bi in range(bs * b):
            cons = cons + [IP_1[bi] + ' + ' + IP_2[bi] + ' >= 1']
            cons = cons + [d1[bi] + ' + ' + IP_2[bi] + ' = 1']
            cons = cons + [d2[bi] + ' + ' + IP_1[bi] + ' = 1']
        return cons

    def genConstraints_forward_round(self, r):
        cons = []
        IP_1 = []
        IP_2 = []
        OP_1 = []
        OP_2 = []
        for pos in self.pos4:
            IP_1.append(Vars_generator.genVars_input_of_round(1, r, pos))
            IP_2.append(Vars_generator.genVars_input_of_round(2, r, pos))
            OP_1.append(Vars_generator.genVars_output_of_round(1, r, pos))
            OP_2.append(Vars_generator.genVars_output_of_round(2, r, pos))
        # - Left and Right Branch
        for i in range(2):
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos3[i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos3[i])
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_1[2 * i]), OSR_1)
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_2[2 * i]), OSR_2)
            # - Separate
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[i])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[i])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[i])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[i])
            MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_Without_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    MC_In_isWhite[bi],
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[i])
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[i])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[i])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[i])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[i])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[i])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[i])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[i])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[i])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[i])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[i])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[i])
            for Coli in range(ColN):
                cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                    column(IMC_SupP_Blue_1, Coli),
                    column(IMC_SupP_Blue_2, Coli),
                    IMC_SupP_Blue_ColExistWhite[Coli],
                    IMC_SupP_Blue_ColAllGray[Coli],
                    column(OMC_SupP_Blue_1, Coli),
                    column(OMC_SupP_Blue_2, Coli),
                    G_IMC_SupP_Blue_SumGray[Coli],
                    G_MC_CD_Blue[Coli]
                )
                cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                    column(IMC_SupP_Red_1, Coli),
                    column(IMC_SupP_Red_2, Coli),
                    IMC_SupP_Red_ColExistWhite[Coli],
                    IMC_SupP_Red_ColAllGray[Coli],
                    column(OMC_SupP_Red_1, Coli),
                    column(OMC_SupP_Red_2, Coli),
                    G_IMC_SupP_Red_SumGray[Coli],
                    G_MC_CD_Red[Coli]
                )
            # - Separate
            IXor_1 = IP_1[2 * i + 1]
            IXor_2 = IP_2[2 * i + 1]
            IXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_Xor(1, r, self.op[2 + i])
            IXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_Xor(2, r, self.op[2 + i])
            IXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_Xor(1, r, self.op[2 + i])
            IXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_Xor(2, r, self.op[2 + i])
            Xor_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[2 + i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_Without_Guess_i(
                    IXor_1[bi],
                    IXor_2[bi],
                    IXor_SupP_Blue_1[bi],
                    IXor_SupP_Blue_2[bi],
                    IXor_SupP_Red_1[bi],
                    IXor_SupP_Red_2[bi],
                    Xor_In_isWhite[bi],
                )
            # - Xor
            OXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_Xor(1, r, self.op[2 + i])
            OXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_Xor(2, r, self.op[2 + i])
            OXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_Xor(1, r, self.op[2 + i])
            OXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_Xor(2, r, self.op[2 + i])
            CD_Xor_Blue = Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, self.op[2 + i])
            CD_Xor_Red = Vars_generator.genVars_Xor_ConsumedDeg_Red(r, self.op[2 + i])
            OXor_isWhite = Vars_generator.genVars_OXor_isWhite(r, self.op[2 + i])
            OXor_SupP_Blue_AND_1 = Vars_generator.genVars_OXor_SupP_Blue_AND(1, r, self.op[2 + i])
            OXor_SupP_Blue_AND_2 = Vars_generator.genVars_OXor_SupP_Blue_AND(2, r, self.op[2 + i])
            OXor_SupP_Blue_OR_1 = Vars_generator.genVars_OXor_SupP_Blue_OR(1, r, self.op[2 + i])
            OXor_SupP_Blue_OR_2 = Vars_generator.genVars_OXor_SupP_Blue_OR(2, r, self.op[2 + i])
            OXor_SupP_Red_AND_1 = Vars_generator.genVars_OXor_SupP_Red_AND(1, r, self.op[2 + i])
            OXor_SupP_Red_AND_2 = Vars_generator.genVars_OXor_SupP_Red_AND(2, r, self.op[2 + i])
            OXor_SupP_Red_OR_1 = Vars_generator.genVars_OXor_SupP_Red_OR(1, r, self.op[2 + i])
            OXor_SupP_Red_OR_2 = Vars_generator.genVars_OXor_SupP_Red_OR(2, r, self.op[2 + i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.genConstrains_of_Xor_i(
                    OMC_SupP_Blue_1[bi],
                    OMC_SupP_Blue_2[bi],
                    OMC_SupP_Red_1[bi],
                    OMC_SupP_Red_2[bi],
                    IXor_SupP_Blue_1[bi],
                    IXor_SupP_Blue_2[bi],
                    IXor_SupP_Red_1[bi],
                    IXor_SupP_Red_2[bi],
                    OXor_SupP_Blue_1[bi],
                    OXor_SupP_Blue_2[bi],
                    OXor_SupP_Red_1[bi],
                    OXor_SupP_Red_2[bi],
                    CD_Xor_Blue[bi],
                    CD_Xor_Red[bi],
                    OXor_isWhite[bi],
                    OXor_SupP_Blue_AND_1[bi],
                    OXor_SupP_Blue_AND_2[bi],
                    OXor_SupP_Blue_OR_1[bi],
                    OXor_SupP_Blue_OR_2[bi],
                    OXor_SupP_Red_AND_1[bi],
                    OXor_SupP_Red_AND_2[bi],
                    OXor_SupP_Red_OR_1[bi],
                    OXor_SupP_Red_OR_2[bi]
                )
            # - Merge
            OXor_1 = Vars_generator.genVars_output_of_Xor(1, r, self.pos3[i])
            OXor_2 = Vars_generator.genVars_output_of_Xor(2, r, self.pos3[i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[bi], OXor_SupP_Red_1[bi]], OXor_1[bi])
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[bi], OXor_SupP_Red_2[bi]], OXor_2[bi])
        # - Middle Branch
        OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos3[2])
        OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos3[2])
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(ShiftRow(IP_1[2])), OSR_1)
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(ShiftRow(IP_2[2])), OSR_2)
        # - Separate
        IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[4])
        IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[4])
        IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[4])
        IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[4])
        MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[4])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Separate_Without_Guess_i(
                OSR_1[bi],
                OSR_2[bi],
                IMC_SupP_Blue_1[bi],
                IMC_SupP_Blue_2[bi],
                IMC_SupP_Red_1[bi],
                IMC_SupP_Red_2[bi],
                MC_In_isWhite[bi],
            )
        # - MixColumn
        IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[4])
        IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[4])
        G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[4])
        G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[4])
        IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[4])
        IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[4])
        G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[4])
        G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[4])
        OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[4])
        OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[4])
        OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[4])
        OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[4])
        for Coli in range(ColN):
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                column(IMC_SupP_Blue_1, Coli),
                column(IMC_SupP_Blue_2, Coli),
                IMC_SupP_Blue_ColExistWhite[Coli],
                IMC_SupP_Blue_ColAllGray[Coli],
                column(OMC_SupP_Blue_1, Coli),
                column(OMC_SupP_Blue_2, Coli),
                G_IMC_SupP_Blue_SumGray[Coli],
                G_MC_CD_Blue[Coli]
            )
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                column(IMC_SupP_Red_1, Coli),
                column(IMC_SupP_Red_2, Coli),
                IMC_SupP_Red_ColExistWhite[Coli],
                IMC_SupP_Red_ColAllGray[Coli],
                column(OMC_SupP_Red_1, Coli),
                column(OMC_SupP_Red_2, Coli),
                G_IMC_SupP_Red_SumGray[Coli],
                G_MC_CD_Red[Coli]
            )
        # - Merge
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]], OP_1[1][bi])
            cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]], OP_2[1][bi])
        # - Link the input and output
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_1[0]), OP_1[3])
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_2[0]), OP_2[3])
        OXor_1 = Vars_generator.genVars_output_of_Xor(1, r, self.pos3[0])
        OXor_2 = Vars_generator.genVars_output_of_Xor(2, r, self.pos3[0])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_1, OP_1[0])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_2, OP_2[0])
        OXor_1 = Vars_generator.genVars_output_of_Xor(1, r, self.pos3[1])
        OXor_2 = Vars_generator.genVars_output_of_Xor(2, r, self.pos3[1])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_1, OP_1[2])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_2, OP_2[2])
        return cons

    def genConstraints_backward_round(self, r):
        cons = []
        IP_1 = []
        IP_2 = []
        OP_1 = []
        OP_2 = []
        for pos in self.pos4:
            IP_1.append(Vars_generator.genVars_input_of_round(1, r, pos))
            IP_2.append(Vars_generator.genVars_input_of_round(2, r, pos))
            OP_1.append(Vars_generator.genVars_output_of_round(1, r, pos))
            OP_2.append(Vars_generator.genVars_output_of_round(2, r, pos))
        # - Left and Right Branch
        for i in range(2):
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos3[i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos3[i])
            if i == 0:
                cons = cons + MITMPreConstraints.equalConstraints(OP_1[3], OSR_1)
                cons = cons + MITMPreConstraints.equalConstraints(OP_2[3], OSR_2)
            if i == 1:
                OSR_1_Middle = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos3[2])
                OSR_2_Middle = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos3[2])
                cons = cons + MITMPreConstraints.equalConstraints(ShiftRow_Inv(OSR_1_Middle), OSR_1)
                cons = cons + MITMPreConstraints.equalConstraints(ShiftRow_Inv(OSR_2_Middle), OSR_2)
            # - Separate
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[i])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[i])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[i])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[i])
            MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_Without_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    MC_In_isWhite[bi],
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[i])
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[i])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[i])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[i])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[i])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[i])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[i])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[i])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[i])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[i])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[i])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[i])
            for Coli in range(ColN):
                cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                    column(IMC_SupP_Blue_1, Coli),
                    column(IMC_SupP_Blue_2, Coli),
                    IMC_SupP_Blue_ColExistWhite[Coli],
                    IMC_SupP_Blue_ColAllGray[Coli],
                    column(OMC_SupP_Blue_1, Coli),
                    column(OMC_SupP_Blue_2, Coli),
                    G_IMC_SupP_Blue_SumGray[Coli],
                    G_MC_CD_Blue[Coli]
                )
                cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                    column(IMC_SupP_Red_1, Coli),
                    column(IMC_SupP_Red_2, Coli),
                    IMC_SupP_Red_ColExistWhite[Coli],
                    IMC_SupP_Red_ColAllGray[Coli],
                    column(OMC_SupP_Red_1, Coli),
                    column(OMC_SupP_Red_2, Coli),
                    G_IMC_SupP_Red_SumGray[Coli],
                    G_MC_CD_Red[Coli]
                )
            # - Separate
            IXor_1 = OP_1[2 * i]
            IXor_2 = OP_2[2 * i]
            IXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_Xor(1, r, self.op[2 + i])
            IXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_Xor(2, r, self.op[2 + i])
            IXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_Xor(1, r, self.op[2 + i])
            IXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_Xor(2, r, self.op[2 + i])
            Xor_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[2 + i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_Without_Guess_i(
                    IXor_1[bi],
                    IXor_2[bi],
                    IXor_SupP_Blue_1[bi],
                    IXor_SupP_Blue_2[bi],
                    IXor_SupP_Red_1[bi],
                    IXor_SupP_Red_2[bi],
                    Xor_In_isWhite[bi],
                )
            # - Xor
            OXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_Xor(1, r, self.op[2 + i])
            OXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_Xor(2, r, self.op[2 + i])
            OXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_Xor(1, r, self.op[2 + i])
            OXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_Xor(2, r, self.op[2 + i])
            CD_Xor_Blue = Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, self.op[2 + i])
            CD_Xor_Red = Vars_generator.genVars_Xor_ConsumedDeg_Red(r, self.op[2 + i])
            OXor_isWhite = Vars_generator.genVars_OXor_isWhite(r, self.op[2 + i])
            OXor_SupP_Blue_AND_1 = Vars_generator.genVars_OXor_SupP_Blue_AND(1, r, self.op[2 + i])
            OXor_SupP_Blue_AND_2 = Vars_generator.genVars_OXor_SupP_Blue_AND(2, r, self.op[2 + i])
            OXor_SupP_Blue_OR_1 = Vars_generator.genVars_OXor_SupP_Blue_OR(1, r, self.op[2 + i])
            OXor_SupP_Blue_OR_2 = Vars_generator.genVars_OXor_SupP_Blue_OR(2, r, self.op[2 + i])
            OXor_SupP_Red_AND_1 = Vars_generator.genVars_OXor_SupP_Red_AND(1, r, self.op[2 + i])
            OXor_SupP_Red_AND_2 = Vars_generator.genVars_OXor_SupP_Red_AND(2, r, self.op[2 + i])
            OXor_SupP_Red_OR_1 = Vars_generator.genVars_OXor_SupP_Red_OR(1, r, self.op[2 + i])
            OXor_SupP_Red_OR_2 = Vars_generator.genVars_OXor_SupP_Red_OR(2, r, self.op[2 + i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.genConstrains_of_Xor_i(
                    OMC_SupP_Blue_1[bi],
                    OMC_SupP_Blue_2[bi],
                    OMC_SupP_Red_1[bi],
                    OMC_SupP_Red_2[bi],
                    IXor_SupP_Blue_1[bi],
                    IXor_SupP_Blue_2[bi],
                    IXor_SupP_Red_1[bi],
                    IXor_SupP_Red_2[bi],
                    OXor_SupP_Blue_1[bi],
                    OXor_SupP_Blue_2[bi],
                    OXor_SupP_Red_1[bi],
                    OXor_SupP_Red_2[bi],
                    CD_Xor_Blue[bi],
                    CD_Xor_Red[bi],
                    OXor_isWhite[bi],
                    OXor_SupP_Blue_AND_1[bi],
                    OXor_SupP_Blue_AND_2[bi],
                    OXor_SupP_Blue_OR_1[bi],
                    OXor_SupP_Blue_OR_2[bi],
                    OXor_SupP_Red_AND_1[bi],
                    OXor_SupP_Red_AND_2[bi],
                    OXor_SupP_Red_OR_1[bi],
                    OXor_SupP_Red_OR_2[bi]
                )
            # - Merge
            OXor_1 = Vars_generator.genVars_output_of_Xor(1, r, self.pos3[i])
            OXor_2 = Vars_generator.genVars_output_of_Xor(2, r, self.pos3[i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[bi], OXor_SupP_Red_1[bi]], OXor_1[bi])
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[bi], OXor_SupP_Red_2[bi]], OXor_2[bi])
        # - Middle Branch
        # - Separate
        IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[4])
        IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[4])
        IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[4])
        IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[4])
        MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[4])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Separate_Without_Guess_i(
                OP_1[1][bi],
                OP_2[1][bi],
                IMC_SupP_Blue_1[bi],
                IMC_SupP_Blue_2[bi],
                IMC_SupP_Red_1[bi],
                IMC_SupP_Red_2[bi],
                MC_In_isWhite[bi],
            )
        # - MixColumn
        IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[4])
        IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[4])
        G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[4])
        G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[4])
        IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[4])
        IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[4])
        G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[4])
        G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[4])
        OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[4])
        OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[4])
        OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[4])
        OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[4])
        for Coli in range(ColN):
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                column(IMC_SupP_Blue_1, Coli),
                column(IMC_SupP_Blue_2, Coli),
                IMC_SupP_Blue_ColExistWhite[Coli],
                IMC_SupP_Blue_ColAllGray[Coli],
                column(OMC_SupP_Blue_1, Coli),
                column(OMC_SupP_Blue_2, Coli),
                G_IMC_SupP_Blue_SumGray[Coli],
                G_MC_CD_Blue[Coli]
            )
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                column(IMC_SupP_Red_1, Coli),
                column(IMC_SupP_Red_2, Coli),
                IMC_SupP_Red_ColExistWhite[Coli],
                IMC_SupP_Red_ColAllGray[Coli],
                column(OMC_SupP_Red_1, Coli),
                column(OMC_SupP_Red_2, Coli),
                G_IMC_SupP_Red_SumGray[Coli],
                G_MC_CD_Red[Coli]
            )
        # - Merge
        OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos3[2])
        OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos3[2])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]], OSR_1[bi])
            cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]], OSR_2[bi])
        # - Link the input and output
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(ShiftRow(IP_1[2])), OSR_1)
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(ShiftRow(IP_2[2])), OSR_2)
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow_Inv(OP_1[3]), IP_1[0])
        cons = cons + MITMPreConstraints.equalConstraints(ShiftRow_Inv(OP_2[3]), IP_2[0])
        OXor_1 = Vars_generator.genVars_output_of_Xor(1, r, self.pos3[0])
        OXor_2 = Vars_generator.genVars_output_of_Xor(2, r, self.pos3[0])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_1, IP_1[1])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_2, IP_2[1])
        OXor_1 = Vars_generator.genVars_output_of_Xor(1, r, self.pos3[1])
        OXor_2 = Vars_generator.genVars_output_of_Xor(2, r, self.pos3[1])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_1, IP_1[3])
        cons = cons + MITMPreConstraints.equalConstraints(OXor_2, IP_2[3])

        return cons

    def genConstraints_link_round(self):
        cons = []
        for r in range(self.TR):
            next_r = (r + 1) % self.TR
            OP_1 = []
            OP_2 = []
            IP_1_next_r = []
            IP_2_next_r = []
            for pos in self.pos4:
                OP_1.append(Vars_generator.genVars_output_of_round(1, r, pos))
                OP_2.append(Vars_generator.genVars_output_of_round(2, r, pos))
            for pos in self.pos4:
                IP_1_next_r.append(Vars_generator.genVars_input_of_round(1, next_r, pos))
                IP_2_next_r.append(Vars_generator.genVars_input_of_round(2, next_r, pos))

            if r == self.TR - 1:
                if self.ini_r <= self.mat_r:
                    # - The last round is backward
                    for i in range(b):
                        if self.output_index[i] == 0:
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 0))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 0))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 1))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 1))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 2), column(IP_1_next_r[0], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 2), column(IP_2_next_r[0], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 3), column(IP_1_next_r[0], 3))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 3), column(IP_2_next_r[0], 3))
                        if self.output_index[i] == 1:
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 0))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 0))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 1))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 1))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 2), column(IP_1_next_r[1], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 2), column(IP_2_next_r[1], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 3), column(IP_1_next_r[1], 3))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 3), column(IP_2_next_r[1], 3))
                        if self.output_index[i] == 2:
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 2))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 2))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 3))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 3))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 0), column(IP_1_next_r[2], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 0), column(IP_2_next_r[2], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 1), column(IP_1_next_r[2], 1))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 1), column(IP_2_next_r[2], 1))
                        if self.output_index[i] == 3:
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 2))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 2))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_1[i], 3))
                            cons = cons + MITMPreConstraints.equalZeroConstraints(column(OP_2[i], 3))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 0), column(IP_1_next_r[3], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 0), column(IP_2_next_r[3], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 1), column(IP_1_next_r[3], 1))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 1), column(IP_2_next_r[3], 1))

                if self.ini_r > self.mat_r:
                    # - The last round is forward
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[0], 0))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[0], 0))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[0], 1))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[0], 1))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[1], 0))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[1], 0))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[1], 1))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[1], 1))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[2], 2))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[2], 2))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[2], 3))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[2], 3))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[3], 2))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[3], 2))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_1_next_r[3], 3))
                    cons = cons + MITMPreConstraints.equalZeroConstraints(column(IP_2_next_r[3], 3))
                    for i in range(b):
                        if self.output_index[i] == 0:
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 2), column(IP_1_next_r[0], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 2), column(IP_2_next_r[0], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 3), column(IP_1_next_r[0], 3))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 3), column(IP_2_next_r[0], 3))
                        if self.output_index[i] == 1:
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 2), column(IP_1_next_r[1], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 2), column(IP_2_next_r[1], 2))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 3), column(IP_1_next_r[1], 3))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 3), column(IP_2_next_r[1], 3))
                        if self.output_index[i] == 2:
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 0), column(IP_1_next_r[2], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 0), column(IP_2_next_r[2], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 1), column(IP_1_next_r[2], 1))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 1), column(IP_2_next_r[2], 1))
                        if self.output_index[i] == 3:
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 0), column(IP_1_next_r[3], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 0), column(IP_2_next_r[3], 0))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_1[i], 1), column(IP_1_next_r[3], 1))
                            cons = cons + MITMPreConstraints.equalConstraints(column(OP_2[i], 1), column(IP_2_next_r[3], 1))
            else:
                for i in range(b):
                    cons = cons + MITMPreConstraints.equalConstraints(OP_1[i], IP_1_next_r[i])
                    cons = cons + MITMPreConstraints.equalConstraints(OP_2[i], IP_2_next_r[i])
        return cons

    def genConstraints_Match(self):
        cons = []
        IP_1 = []
        IP_2 = []
        OP_1 = []
        OP_2 = []
        for pos in self.pos4:
            IP_1.append(Vars_generator.genVars_input_of_round(1, self.mat_r, pos))
            IP_2.append(Vars_generator.genVars_input_of_round(2, self.mat_r, pos))
            OP_1.append(Vars_generator.genVars_output_of_round(1, self.mat_r, pos))
            OP_2.append(Vars_generator.genVars_output_of_round(2, self.mat_r, pos))
        # - Left Branch
        G_Match_Counter = Vars_generator.genVar_Match_Counted(self.pos3[0])
        OSR_1 = ShiftRow(IP_1[0])
        OSR_2 = ShiftRow(IP_2[0])
        Match_IMC_isWhite = Vars_generator.genVars_Match_IMC_isWhite(self.pos3[0])
        Match_IP_isWhite = Vars_generator.genVars_Match_IP_isWhite(self.pos3[0])
        Match_OP_isWhite = Vars_generator.genVars_Match_OP_isWhite(self.pos3[0])
        Match_OXor_isWhite = Vars_generator.genVars_Match_OXor_isWhite(self.pos3[0])
        Match_Exist = Vars_generator.genVars_Match_Exist(self.pos3[0])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OSR_1[bi], OSR_2[bi]], Match_IMC_isWhite[bi])
            cons = cons + MITMPreConstraints.Determine_Allzero([IP_1[1][bi], IP_2[1][bi]], Match_IP_isWhite[bi])
            cons = cons + MITMPreConstraints.Determine_Allzero([OP_1[0][bi], OP_2[0][bi]], Match_OP_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_ExistOne([Match_IP_isWhite[bi], Match_OP_isWhite[bi]], Match_OXor_isWhite[bi])
        for Coli in range(ColN):
            Match_IMC_isWhite_Coli = column(Match_IMC_isWhite, Coli)
            Match_OXor_isWhite_Coli = column(Match_OXor_isWhite, Coli)
            cons = cons + [str(RowN + 1) + ' ' + Match_Exist[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OXor_isWhite_Coli) + ' <= ' + str(SumIOMC)]
            cons = cons + [str(RowN) + ' ' + Match_Exist[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OXor_isWhite_Coli) + ' >= ' + str(RowN)]
            cons = cons + [Match_Exist[Coli] + ' = 1 -> ' + G_Match_Counter[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OXor_isWhite_Coli) + ' = ' + str(RowN)]
            cons = cons + [Match_Exist[Coli] + ' = 0 -> ' + G_Match_Counter[Coli] + ' = 0']

        # - Right Branch
        G_Match_Counter = Vars_generator.genVar_Match_Counted(self.pos3[1])
        OSR_1 = ShiftRow(IP_1[2])
        OSR_2 = ShiftRow(IP_2[2])
        Match_IMC_isWhite = Vars_generator.genVars_Match_IMC_isWhite(self.pos3[1])
        Match_IP_isWhite = Vars_generator.genVars_Match_IP_isWhite(self.pos3[1])
        Match_OP_isWhite = Vars_generator.genVars_Match_OP_isWhite(self.pos3[1])
        Match_OXor_isWhite = Vars_generator.genVars_Match_OXor_isWhite(self.pos3[1])
        Match_Exist = Vars_generator.genVars_Match_Exist(self.pos3[1])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OSR_1[bi], OSR_2[bi]], Match_IMC_isWhite[bi])
            cons = cons + MITMPreConstraints.Determine_Allzero([IP_1[3][bi], IP_2[3][bi]], Match_IP_isWhite[bi])
            cons = cons + MITMPreConstraints.Determine_Allzero([OP_1[2][bi], OP_2[2][bi]], Match_OP_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_ExistOne([Match_IP_isWhite[bi], Match_OP_isWhite[bi]], Match_OXor_isWhite[bi])
        for Coli in range(ColN):
            Match_IMC_isWhite_Coli = column(Match_IMC_isWhite, Coli)
            Match_OXor_isWhite_Coli = column(Match_OXor_isWhite, Coli)
            cons = cons + [str(RowN + 1) + ' ' + Match_Exist[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OXor_isWhite_Coli) + ' <= ' + str(SumIOMC)]
            cons = cons + [str(RowN) + ' ' + Match_Exist[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OXor_isWhite_Coli) + ' >= ' + str(RowN)]
            cons = cons + [Match_Exist[Coli] + ' = 1 -> ' + G_Match_Counter[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OXor_isWhite_Coli) + ' = ' + str(RowN)]
            cons = cons + [Match_Exist[Coli] + ' = 0 -> ' + G_Match_Counter[Coli] + ' = 0']

        # - Middle Branch
        G_Match_Counter = Vars_generator.genVar_Match_Counted(self.pos3[2])
        OSR_1 = ShiftRow(ShiftRow(IP_1[2]))
        OSR_2 = ShiftRow(ShiftRow(IP_2[2]))
        Match_IMC_isWhite = Vars_generator.genVars_Match_IMC_isWhite(self.pos3[2])
        Match_OP_isWhite = Vars_generator.genVars_Match_OP_isWhite(self.pos3[2])
        Match_Exist = Vars_generator.genVars_Match_Exist(self.pos3[2])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OSR_1[bi], OSR_2[bi]], Match_IMC_isWhite[bi])
            cons = cons + MITMPreConstraints.Determine_Allzero([OP_1[1][bi], OP_2[1][bi]], Match_OP_isWhite[bi])
        for Coli in range(ColN):
            Match_IMC_isWhite_Coli = column(Match_IMC_isWhite, Coli)
            Match_OP_isWhite_Coli = column(Match_OP_isWhite, Coli)
            cons = cons + [str(RowN + 1) + ' ' + Match_Exist[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OP_isWhite_Coli) + ' <= ' + str(SumIOMC)]
            cons = cons + [str(RowN) + ' ' + Match_Exist[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OP_isWhite_Coli) + ' >= ' + str(RowN)]
            cons = cons + [Match_Exist[Coli] + ' = 1 -> ' + G_Match_Counter[Coli] + ' + ' + BasicTools.plusTerms(Match_IMC_isWhite_Coli + Match_OP_isWhite_Coli) + ' = ' + str(RowN)]
            cons = cons + [Match_Exist[Coli] + ' = 0 -> ' + G_Match_Counter[Coli] + ' = 0']

        # - Direct Match
        G_Match_Counter = Vars_generator.genVars_Match_Counted_direct()
        OSR_1 = ShiftRow(IP_1[0])
        OSR_2 = ShiftRow(IP_2[0])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Match_direct([OSR_1[bi], OSR_2[bi]], [OP_1[3][bi], OP_2[3][bi]], G_Match_Counter[bi])
        return cons

    def genConstraints_additional(self):
        cons = []
        CD_Blue = []
        CD_Red = []
        for r in range(self.TR):
            if r != self.mat_r:
                CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[0])
                CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[1])
                CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[4])
                CD_Blue = CD_Blue + Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, self.op[2])
                CD_Blue = CD_Blue + Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, self.op[3])
                CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[0])
                CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[1])
                CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[4])
                CD_Red = CD_Red + Vars_generator.genVars_Xor_ConsumedDeg_Red(r, self.op[2])
                CD_Red = CD_Red + Vars_generator.genVars_Xor_ConsumedDeg_Red(r, self.op[3])

        d1 = Vars_generator.genVars_degree__forward()
        d2 = Vars_generator.genVars_degree_backward()

        Deg1 = 'GDeg1'
        Deg2 = 'GDeg2'

        if len(CD_Blue) > 0:
            cons = cons + [Deg1 + ' - ' + BasicTools.minusTerms(d1) + ' + ' + BasicTools.plusTerms(CD_Blue) + ' = 0']
        else:
            cons = cons + [Deg1 + ' - ' + BasicTools.minusTerms(d1) + ' = 0']
        if len(CD_Red) > 0:
            cons = cons + [Deg2 + ' - ' + BasicTools.minusTerms(d2) + ' + ' + BasicTools.plusTerms(CD_Red) + ' = 0']
        else:
            cons = cons + [Deg2 + ' - ' + BasicTools.minusTerms(d2) + ' = 0']

        cons = cons + [Deg1 + ' >= 1']
        cons = cons + [Deg2 + ' >= 1']

        # DoF = 'GDOF'
        # DoB = 'GDOB'
        # cons = cons + [DoF + ' + ' + BasicTools.plusTerms(d1) + ' = ' + str(h)]
        # cons = cons + [DoB + ' + ' + BasicTools.plusTerms(d2) + ' = ' + str(h)]
        cons = cons + [BasicTools.plusTerms(d1) + ' <= ' + str(h - 1)]
        cons = cons + [BasicTools.plusTerms(d2) + ' <= ' + str(h - 1)]
        # cons = cons + [DoF + ' >= 1']
        # cons = cons + [DoB + ' >= 1']

        GM = 'GMat'
        G_Match_counter = []
        G_Match_counter = G_Match_counter + Vars_generator.genVar_Match_Counted(self.pos3[0])
        G_Match_counter = G_Match_counter + Vars_generator.genVar_Match_Counted(self.pos3[1])
        G_Match_counter = G_Match_counter + Vars_generator.genVar_Match_Counted(self.pos3[2])
        G_Match_counter = G_Match_counter + Vars_generator.genVars_Match_Counted_direct()
        cons = cons + [GM + ' - ' + BasicTools.minusTerms(G_Match_counter) + ' = 0']
        cons = cons + [GM + ' >= 1']
        return cons

    def genConstraints_total(self):
        cons = []
        cons = cons + self.genConstraints_initial_degree()
        if self.ini_r <= self.mat_r:
            for r in range(self.ini_r, self.mat_r):
                cons = cons + self.genConstraints_forward_round(r)
            for r in range(0, self.ini_r):
                cons = cons + self.genConstraints_backward_round(r)
            for r in range(self.mat_r + 1, self.TR):
                cons = cons + self.genConstraints_backward_round(r)
        if self.ini_r > self.mat_r:
            for r in range(self.ini_r, self.TR):
                cons = cons + self.genConstraints_forward_round(r)
            for r in range(0, self.mat_r):
                cons = cons + self.genConstraints_forward_round(r)
            for r in range(self.mat_r + 1, self.ini_r):
                cons = cons + self.genConstraints_backward_round(r)
        cons = cons + self.genConstraints_Match()
        cons = cons + self.genConstraints_link_round()
        cons = cons + self.genConstraints_additional()
        return cons

    def genModel(self, filename):
        V = set([])
        cons = []
        cons = cons + self.genConstraints_total()

        # cons = cons + ['GDeg1 + GDeg2 >= 32']
        cons = cons + ['GObj - GDeg1 <= 0']
        cons = cons + ['GObj - GDeg2 <= 0']
        cons = cons + ['GObj - GMat <= 0']
        # cons = cons + ['GObj - GDOF <= 0']
        # cons = cons + ['GObj - GDOB <= 0']

        V = BasicTools.getVariables_From_Constraints(cons)

        with open(filename + ".lp", "w") as fid:
            fid.write('Maximize' + '\n')
            fid.write('GObj' + '\n')
            fid.write('\n')
            fid.write('Subject To')
            fid.write('\n')
            for c in cons:
                fid.write(c)
                fid.write('\n')

            GV = []
            BV = []
            for v in V:
                if v[0] == 'G':
                    GV.append(v)
                else:
                    BV.append(v)

            fid.write('Binary' + '\n')
            for bv in BV:
                fid.write(bv + '\n')

            fid.write('Generals' + '\n')
            for gv in GV:
                fid.write(gv + '\n')


if __name__ == '__main__':
    TR = int(sys.argv[1])
    root = f'./Model/TR{TR}'
    if not os.path.exists(root):
        os.makedirs(root)
    with open(f"./Model/Result_{TR}.txt", "w") as rd:
        rd.write('TR, ini_r, mat_r: d1, d2, m' + '\n')
        for mat_r in [0, TR - 1, TR - 2, TR - 3]:
            for ini_r in range(TR):
                if ini_r != mat_r:
                    filename = f'./Model/TR{TR}/inir{ini_r}_matr{mat_r}'
                    A = Constraints_generator(TR, ini_r, mat_r)
                    A.genModel(filename)
                    Model = read(filename + '.lp')
                    # Model.setParam('TimeLimit', 180 * 60)
                    Model.optimize()

                    if Model.SolCount == 0:
                        pass
                    else:
                        Model.write(filename + '.sol')
                        solFile = open(filename + '.sol', 'r')
                        Sol = dict()

                        for line in solFile:
                            if line[0] != '#':
                                temp = line
                                # temp = temp.replace('-', ' ')
                                temp = temp.split()
                                Sol[temp[0]] = int(eval(temp[1]))
                        rd.write(str(TR) + ',' + str(ini_r) + ',' + str(mat_r) + ':')
                        rd.write(str(Sol['GDeg1']) + ',' + str(Sol['GDeg2']) + ',' + str(Sol['GMat']) + '\n')
                        rd.flush()
