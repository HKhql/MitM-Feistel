

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

    @staticmethod
    def genVars_Guess_Blue(r, op):
        return [f'P_Guess_Blue_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Guess_Red(r, op):
        return [f'P_Guess_Red_r{r}_{op}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Guess_Both(r, op):
        return [f'P_Guess_Both_r{r}_{op}_{j}' for j in range(bs)]

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

    # Match
    @staticmethod
    def genVars_Match_isWhite(r):
        return [f'Match_isWhite_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_Exist(i):
        return [f'Match_Exist_link{i}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_Counter(i):
        return [f'Match_Counter_{i}_{j}' for j in range(bs)]

    @staticmethod
    def genVasr_number_of_initial_degree():
        return 'Indicator_degree_less'


class Constraints_generator:
    def __init__(self, total_round, initial_round, matching_round):
        self.ini_r = initial_round
        self.mat_r = matching_round
        self.TR = total_round
        self.pos6 = ['LL', 'LR', 'ML', 'MR', 'RL', 'RR']
        self.pos3 = ['L', 'M', 'R']
        self.op_MC = ['MCLL', 'MCLR', 'MCML', 'MCMR', 'MCRL', 'MCRR']
        self.op_Xor = ['XorL', 'XorM', 'XorR']
        self.IO_index = [3, 0, 1, 4, 5, 2]
        self.output_index = output_index()[self.TR]
        self.link = link()
        self.last_round_forward = (self.ini_r > self.mat_r)

    def genConstraints_initial_degree(self):
        cons = []
        IP_1 = []
        IP_2 = []
        d1 = Vars_generator.genVars_degree__forward()
        d2 = Vars_generator.genVars_degree_backward()
        for pos in self.pos6:
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
        for pos in self.pos6:
            IP_1.append(Vars_generator.genVars_input_of_round(1, r, pos))
            IP_2.append(Vars_generator.genVars_input_of_round(2, r, pos))
            OP_1.append(Vars_generator.genVars_output_of_round(1, r, pos))
            OP_2.append(Vars_generator.genVars_output_of_round(2, r, pos))
        # - ShiftRow
        for i in range(3):
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos6[2 * i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos6[2 * i])
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_1[2 * i]), OSR_1)
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_2[2 * i]), OSR_2)
        for i in range(3):
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos6[2 * i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos6[2 * i])
            # - Separate
            In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op_MC[2 * i])
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op_MC[2 * i])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op_MC[2 * i])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op_MC[2 * i])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op_MC[2 * i])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op_MC[2 * i])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op_MC[2 * i])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op_MC[2 * i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op_MC[2 * i])
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op_MC[2 * i])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op_MC[2 * i])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op_MC[2 * i])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op_MC[2 * i])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op_MC[2 * i])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op_MC[2 * i])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op_MC[2 * i])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op_MC[2 * i])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op_MC[2 * i])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op_MC[2 * i])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op_MC[2 * i])
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
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos6[2 * i + 1])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos6[2 * i + 1])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]], ShiftRow_Inv(OSR_1)[bi])
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]], ShiftRow_Inv(OSR_2)[bi])
            # - Separate
            In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op_MC[2 * i + 1])
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op_MC[2 * i + 1])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op_MC[2 * i + 1])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op_MC[2 * i + 1])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op_MC[2 * i + 1])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op_MC[2 * i + 1])
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op_MC[2 * i + 1])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op_MC[2 * i + 1])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op_MC[2 * i + 1])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op_MC[2 * i + 1])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op_MC[2 * i + 1])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op_MC[2 * i + 1])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op_MC[2 * i + 1])
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
            IXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_Xor(1, r, self.op_Xor[i])
            IXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_Xor(2, r, self.op_Xor[i])
            IXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_Xor(1, r, self.op_Xor[i])
            IXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_Xor(2, r, self.op_Xor[i])
            Xor_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op_Xor[i])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op_Xor[i])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op_Xor[i])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op_Xor[i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    IXor_1[bi],
                    IXor_2[bi],
                    IXor_SupP_Blue_1[bi],
                    IXor_SupP_Blue_2[bi],
                    IXor_SupP_Red_1[bi],
                    IXor_SupP_Red_2[bi],
                    Xor_In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - Xor
            OXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_Xor(1, r, self.op_Xor[i])
            OXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_Xor(2, r, self.op_Xor[i])
            OXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_Xor(1, r, self.op_Xor[i])
            OXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_Xor(2, r, self.op_Xor[i])
            CD_Xor_Blue = Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, self.op_Xor[i])
            CD_Xor_Red = Vars_generator.genVars_Xor_ConsumedDeg_Red(r, self.op_Xor[i])
            OXor_isWhite = Vars_generator.genVars_OXor_isWhite(r, self.op_Xor[i])
            OXor_SupP_Blue_AND_1 = Vars_generator.genVars_OXor_SupP_Blue_AND(1, r, self.op_Xor[i])
            OXor_SupP_Blue_AND_2 = Vars_generator.genVars_OXor_SupP_Blue_AND(2, r, self.op_Xor[i])
            OXor_SupP_Blue_OR_1 = Vars_generator.genVars_OXor_SupP_Blue_OR(1, r, self.op_Xor[i])
            OXor_SupP_Blue_OR_2 = Vars_generator.genVars_OXor_SupP_Blue_OR(2, r, self.op_Xor[i])
            OXor_SupP_Red_AND_1 = Vars_generator.genVars_OXor_SupP_Red_AND(1, r, self.op_Xor[i])
            OXor_SupP_Red_AND_2 = Vars_generator.genVars_OXor_SupP_Red_AND(2, r, self.op_Xor[i])
            OXor_SupP_Red_OR_1 = Vars_generator.genVars_OXor_SupP_Red_OR(1, r, self.op_Xor[i])
            OXor_SupP_Red_OR_2 = Vars_generator.genVars_OXor_SupP_Red_OR(2, r, self.op_Xor[i])
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
            OXor_1 = OP_1[self.IO_index[2 * i + 1]]
            OXor_2 = OP_2[self.IO_index[2 * i + 1]]
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[bi], OXor_SupP_Red_1[bi]], OXor_1[bi])
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[bi], OXor_SupP_Red_2[bi]], OXor_2[bi])
        # - Link
        for i in range(3):
            cons = cons + MITMPreConstraints.equalConstraints(IP_1[2 * i], OP_1[self.IO_index[2 * i]])
            cons = cons + MITMPreConstraints.equalConstraints(IP_2[2 * i], OP_2[self.IO_index[2 * i]])
        return cons

    def genConstraints_backward_round(self, r):
        cons = []
        IP_1 = []
        IP_2 = []
        OP_1 = []
        OP_2 = []
        for pos in self.pos6:
            IP_1.append(Vars_generator.genVars_input_of_round(1, r, pos))
            IP_2.append(Vars_generator.genVars_input_of_round(2, r, pos))
            OP_1.append(Vars_generator.genVars_output_of_round(1, r, pos))
            OP_2.append(Vars_generator.genVars_output_of_round(2, r, pos))
        # - ShiftRow
        for i in range(3):
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos6[2 * i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos6[2 * i])
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_1[2 * i]), OSR_1)
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_2[2 * i]), OSR_2)
        for i in range(3):
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos6[2 * i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos6[2 * i])
            # - Separate
            In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op_MC[2 * i])
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op_MC[2 * i])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op_MC[2 * i])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op_MC[2 * i])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op_MC[2 * i])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op_MC[2 * i])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op_MC[2 * i])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op_MC[2 * i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op_MC[2 * i])
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op_MC[2 * i])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op_MC[2 * i])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op_MC[2 * i])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op_MC[2 * i])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op_MC[2 * i])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op_MC[2 * i])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op_MC[2 * i])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op_MC[2 * i])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op_MC[2 * i])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op_MC[2 * i])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op_MC[2 * i])
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
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r, self.pos6[2 * i + 1])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r, self.pos6[2 * i + 1])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]],
                                                                  ShiftRow_Inv(OSR_1)[bi]
                                                                  )
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]],
                                                                  ShiftRow_Inv(OSR_2)[bi]
                                                                  )
            # - Separate
            In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op_MC[2 * i + 1])
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op_MC[2 * i + 1])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op_MC[2 * i + 1])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op_MC[2 * i + 1])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op_MC[2 * i + 1])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r,
                                                                                            self.op_MC[2 * i + 1]
                                                                                            )
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op_MC[2 * i + 1])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op_MC[2 * i + 1])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op_MC[2 * i + 1])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op_MC[2 * i + 1])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op_MC[2 * i + 1])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op_MC[2 * i + 1])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op_MC[2 * i + 1])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op_MC[2 * i + 1])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op_MC[2 * i + 1])
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
            IXor_1 = OP_1[self.IO_index[2 * i + 1]]
            IXor_2 = OP_2[self.IO_index[2 * i + 1]]
            IXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_Xor(1, r, self.op_Xor[i])
            IXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_Xor(2, r, self.op_Xor[i])
            IXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_Xor(1, r, self.op_Xor[i])
            IXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_Xor(2, r, self.op_Xor[i])
            Xor_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op_Xor[i])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op_Xor[i])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op_Xor[i])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op_Xor[i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    IXor_1[bi],
                    IXor_2[bi],
                    IXor_SupP_Blue_1[bi],
                    IXor_SupP_Blue_2[bi],
                    IXor_SupP_Red_1[bi],
                    IXor_SupP_Red_2[bi],
                    Xor_In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - Xor
            OXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_Xor(1, r, self.op_Xor[i])
            OXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_Xor(2, r, self.op_Xor[i])
            OXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_Xor(1, r, self.op_Xor[i])
            OXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_Xor(2, r, self.op_Xor[i])
            CD_Xor_Blue = Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, self.op_Xor[i])
            CD_Xor_Red = Vars_generator.genVars_Xor_ConsumedDeg_Red(r, self.op_Xor[i])
            OXor_isWhite = Vars_generator.genVars_OXor_isWhite(r, self.op_Xor[i])
            OXor_SupP_Blue_AND_1 = Vars_generator.genVars_OXor_SupP_Blue_AND(1, r, self.op_Xor[i])
            OXor_SupP_Blue_AND_2 = Vars_generator.genVars_OXor_SupP_Blue_AND(2, r, self.op_Xor[i])
            OXor_SupP_Blue_OR_1 = Vars_generator.genVars_OXor_SupP_Blue_OR(1, r, self.op_Xor[i])
            OXor_SupP_Blue_OR_2 = Vars_generator.genVars_OXor_SupP_Blue_OR(2, r, self.op_Xor[i])
            OXor_SupP_Red_AND_1 = Vars_generator.genVars_OXor_SupP_Red_AND(1, r, self.op_Xor[i])
            OXor_SupP_Red_AND_2 = Vars_generator.genVars_OXor_SupP_Red_AND(2, r, self.op_Xor[i])
            OXor_SupP_Red_OR_1 = Vars_generator.genVars_OXor_SupP_Red_OR(1, r, self.op_Xor[i])
            OXor_SupP_Red_OR_2 = Vars_generator.genVars_OXor_SupP_Red_OR(2, r, self.op_Xor[i])
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
            OXor_1 = IP_1[2 * i + 1]
            OXor_2 = IP_2[2 * i + 1]
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[bi], OXor_SupP_Red_1[bi]],
                                                                  OXor_1[bi]
                                                                  )
                cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[bi], OXor_SupP_Red_2[bi]],
                                                                  OXor_2[bi]
                                                                  )
        # - Link
        for i in range(3):
            cons = cons + MITMPreConstraints.equalConstraints(IP_1[2 * i], OP_1[self.IO_index[2 * i]])
            cons = cons + MITMPreConstraints.equalConstraints(IP_2[2 * i], OP_2[self.IO_index[2 * i]])
        return cons

    def genConstraints_link_round(self):
        cons = []
        for r in range(self.TR):
            if r != self.mat_r:
                IP_next_r_1 = []
                IP_next_r_2 = []
                OP_1 = []
                OP_2 = []
                for pos in self.pos6:
                    OP_1.append(Vars_generator.genVars_output_of_round(1, r, pos))
                    OP_2.append(Vars_generator.genVars_output_of_round(2, r, pos))
                    IP_next_r_1.append(Vars_generator.genVars_input_of_round(1, (r + 1) % self.TR, pos))
                    IP_next_r_2.append(Vars_generator.genVars_input_of_round(2, (r + 1) % self.TR, pos))
                if r != self.TR - 1:
                    for i in range(6):
                        cons = cons + MITMPreConstraints.equalConstraints(OP_1[i], IP_next_r_1[i])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_2[i], IP_next_r_2[i])
                else:
                    if self.last_round_forward:
                        for i in range(2, 6):
                            cons = cons + MITMPreConstraints.equalZeroConstraints(IP_next_r_1[i])
                            cons = cons + MITMPreConstraints.equalZeroConstraints(IP_next_r_2[i])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_1[self.output_index[0]], IP_next_r_1[0])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_2[self.output_index[0]], IP_next_r_2[0])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_1[self.output_index[1]], IP_next_r_1[1])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_2[self.output_index[1]], IP_next_r_2[1])
                    else:
                        for i in range(6):
                            if i not in self.output_index:
                                cons = cons + MITMPreConstraints.equalZeroConstraints(OP_1[i])
                                cons = cons + MITMPreConstraints.equalZeroConstraints(OP_1[i])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_1[self.output_index[0]], IP_next_r_1[0])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_2[self.output_index[0]], IP_next_r_2[0])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_1[self.output_index[1]], IP_next_r_1[1])
                        cons = cons + MITMPreConstraints.equalConstraints(OP_2[self.output_index[1]], IP_next_r_2[1])
        return cons

    def genConstraints_Match_round(self):
        cons = []
        for i in range(3):
            # - ShiftRow
            IP_1 = Vars_generator.genVars_input_of_round(1, self.mat_r, self.pos6[2 * i])
            IP_2 = Vars_generator.genVars_input_of_round(2, self.mat_r, self.pos6[2 * i])
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, self.mat_r, self.pos6[2 * i])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, self.mat_r, self.pos6[2 * i])
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_1), OSR_1)
            cons = cons + MITMPreConstraints.equalConstraints(ShiftRow(IP_2), OSR_2)
            # - Separate
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, self.mat_r, self.op_MC[2 * i])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, self.mat_r, self.op_MC[2 * i])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, self.mat_r, self.op_MC[2 * i])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, self.mat_r, self.op_MC[2 * i])
            MC_In_isWhite = Vars_generator.genVars_In_isWhite(self.mat_r, self.op_MC[2 * i])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(self.mat_r, self.op_MC[2 * i])
            Guess_Red = Vars_generator.genVars_Guess_Red(self.mat_r, self.op_MC[2 * i])
            Guess_Both = Vars_generator.genVars_Guess_Both(self.mat_r, self.op_MC[2 * i])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    OSR_1[bi],
                    OSR_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    MC_In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # - MixColumn
            IMC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(self.mat_r, self.op_MC[2 * i])
            IMC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(self.mat_r, self.op_MC[2 * i])
            G_IMC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(self.mat_r, self.op_MC[2 * i])
            G_MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(self.mat_r, self.op_MC[2 * i])
            IMC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(self.mat_r, self.op_MC[2 * i])
            IMC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(self.mat_r, self.op_MC[2 * i])
            G_IMC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(self.mat_r, self.op_MC[2 * i])
            G_MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(self.mat_r, self.op_MC[2 * i])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, self.mat_r, self.op_MC[2 * i])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, self.mat_r, self.op_MC[2 * i])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, self.mat_r, self.op_MC[2 * i])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, self.mat_r, self.op_MC[2 * i])
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
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, self.mat_r, self.pos6[2 * i + 1])
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, self.mat_r, self.pos6[2 * i + 1])
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]], ShiftRow_Inv(OSR_1)[bi])
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]], ShiftRow_Inv(OSR_2)[bi])
        return cons

    def genConstraints_Match(self):
        cons = []
        for i in range(2):
            for r_op in self.link[i]:
                Match_isWhite = Vars_generator.genVars_Match_isWhite(r_op[0])
                OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r_op[0], r_op[1])
                OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r_op[0], r_op[1])
                for bi in range(bs):
                    cons = cons + MITMPreConstraints.Determine_Allzero([OSR_1[bi], OSR_2[bi]], Match_isWhite[bi])

        for i in range(2):
            Match_Exist = Vars_generator.genVars_Match_Exist(i)
            for bi in range(bs):
                Match_link = []
                for r_op in self.link[i]:
                    Match_link.append(Vars_generator.genVars_Match_isWhite(r_op[0])[bi])
                cons = cons + BasicTools.OR_(Match_link, Match_Exist[bi])
        for i in range(2):
            G_Match = Vars_generator.genVars_Match_Counter(i)
            Match_Exist = Vars_generator.genVars_Match_Exist(i)
            for bi in range(bs):
                # X = []
                # Y = []
                # for r_op in self.link[i]:
                #     X.append(Vars_generator.genVars_output_of_ShiftRow(1, r_op[0], r_op[1])[bi])
                #     Y.append(Vars_generator.genVars_output_of_ShiftRow(2, r_op[0], r_op[1])[bi])
                cons = cons + [Match_Exist[bi] + ' = 0 -> ' + G_Match[bi] + ' = 1']
                # pre = f'{Match_Exist[bi]} = 1 -> '
                cons = cons + [Match_Exist[bi] + ' = 1 -> ' + G_Match[bi] + ' = 0']
        return cons

    def genConstraints_additional(self):
        cons = []
        CD_Blue = []
        CD_Red = []
        Guess_Blue = []
        Guess_Red = []
        Guess_Both = []
        for r in range(self.TR):
            if r != self.mat_r:
                for op in self.op_MC:
                    CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, op)
                    CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, op)
                    Guess_Blue = Guess_Blue + Vars_generator.genVars_Guess_Blue(r, op)
                    Guess_Red = Guess_Red + Vars_generator.genVars_Guess_Red(r, op)
                    Guess_Both = Guess_Both + Vars_generator.genVars_Guess_Both(r, op)
                for op in self.op_Xor:
                    CD_Blue = CD_Blue + Vars_generator.genVars_Xor_ConsumedDeg_Blue(r, op)
                    CD_Red = CD_Red + Vars_generator.genVars_Xor_ConsumedDeg_Red(r, op)
                    Guess_Blue = Guess_Blue + Vars_generator.genVars_Guess_Blue(r, op)
                    Guess_Red = Guess_Red + Vars_generator.genVars_Guess_Red(r, op)
                    Guess_Both = Guess_Both + Vars_generator.genVars_Guess_Both(r, op)
            else:
                for i in range(3):
                    CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op_MC[2 * i])
                    CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op_MC[2 * i])
                    Guess_Blue = Guess_Blue + Vars_generator.genVars_Guess_Blue(r, self.op_MC[2 * i])
                    Guess_Red = Guess_Red + Vars_generator.genVars_Guess_Red(r, self.op_MC[2 * i])
                    Guess_Both = Guess_Both + Vars_generator.genVars_Guess_Both(r, self.op_MC[2 * i])

        d1 = Vars_generator.genVars_degree__forward()
        d2 = Vars_generator.genVars_degree_backward()

        Deg1 = 'GDeg1'
        Deg2 = 'GDeg2'

        if len(CD_Blue + Guess_Red) > 0:
            cons = cons + [Deg1 + ' - ' + BasicTools.minusTerms(d1) + ' + ' + BasicTools.plusTerms(CD_Blue + Guess_Red) + ' = 0']
        else:
            cons = cons + [Deg1 + ' - ' + BasicTools.minusTerms(d1) + ' = 0']
        if len(CD_Red + Guess_Blue) > 0:
            cons = cons + [Deg2 + ' - ' + BasicTools.minusTerms(d2) + ' + ' + BasicTools.plusTerms(CD_Red + Guess_Blue) + ' = 0']
        else:
            cons = cons + [Deg2 + ' - ' + BasicTools.minusTerms(d2) + ' = 0']

        cons = cons + [Deg1 + ' >= 1']
        cons = cons + [Deg2 + ' >= 1']

        GM = 'GMat'
        G_Match_counter = []
        G_Match_counter = G_Match_counter + Vars_generator.genVars_Match_Counter(0)
        G_Match_counter = G_Match_counter + Vars_generator.genVars_Match_Counter(1)
        cons = cons + [GM + ' - ' + BasicTools.minusTerms(G_Match_counter) + ' + ' + BasicTools.plusTerms(Guess_Blue + Guess_Red + Guess_Both) + ' = 0']
        cons = cons + [GM + ' >= 1']

        Indicator_degree = Vars_generator.genVasr_number_of_initial_degree()
        cons = cons + [str(h + 1) + ' ' + Indicator_degree + ' + ' + BasicTools.plusTerms(d1 + d2) + ' <= ' + str(2 * h)]
        cons = cons + [str(h) + ' ' + Indicator_degree + ' + ' + BasicTools.plusTerms(d1 + d2) + ' >= ' + str(h)]
        cons = cons + [Indicator_degree + ' = 1 -> ' +
                       BasicTools.plusTerms(Guess_Blue) + ' - ' + BasicTools.minusTerms(d2) + ' <= -1']
        cons = cons + [Indicator_degree + ' = 1 -> ' +
                       BasicTools.plusTerms(Guess_Red) + ' - ' + BasicTools.minusTerms(d1) + ' <= -1']
        cons = cons + [Indicator_degree + ' = 0 -> ' +
                       BasicTools.plusTerms(d1 + Guess_Blue + Guess_Both) + ' <= ' + str(h - 1)]
        cons = cons + [Indicator_degree + ' = 0 -> ' +
                       BasicTools.plusTerms(d2 + Guess_Red + Guess_Both) + ' <= ' + str(h - 1)]
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
        cons = cons + self.genConstraints_Match_round()
        cons = cons + self.genConstraints_link_round()
        cons = cons + self.genConstraints_Match()
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
    # TR = 12
    TR = int(sys.argv[1])
    root = f'./Model/TR{TR}'
    if not os.path.exists(root):
        os.makedirs(root)
    with open(f"./Model/Result_{TR}.txt", "w") as rd:
        rd.write('TR, ini_r, mat_r: d1, d2, m' + '\n')
        for mat_r in [0, 1, 2, TR - 1, TR - 2, TR - 3]:
            for ini_r in range(TR):
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























