

import sys
from utils import *
from gurobipy import *
import os


class Vars_generator:
    @staticmethod
    def genVars_input_of_round(i, r, pos):
        return [f'IP_{i}_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_output_of_round(i, r, pos):
        return [f'OP_{i}_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_output_of_ShiftRow(i, r):
        return [f'OSR_{i}_r{r}_{j}' for j in range(bs)]

    # @staticmethod
    # def genVars_output_of_Xor(i, r):
    #     return [f'OXor_{i}_r{r}_{j}' for j in range(bs)]

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
    def genVars_Xor_ConsumedDeg_Blue(r):
        return [f'CD_Xor_Blue_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Xor_ConsumedDeg_Red(r):
        return [f'CD_Xor_Red__r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Blue_input_of_Xor(i, r):
        return [f'IXor_SupP_Blue_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Red_input_of_Xor(i, r):
        return [f'IXor_SupP_Red_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Blue_output_of_Xor(i, r):
        return [f'OXor_SupP_Blue_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_SupP_Red_output_of_Xor(i, r):
        return [f'OXor_SupP_Red_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_isWhite(r):
        return [f'OXor_isWhite_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Blue_AND(i, r):
        return [f'OXor_SupP_Blue_AND_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Blue_OR(i, r):
        return [f'OXor_SupP_Blue_OR_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Red_AND(i, r):
        return [f'OXor_SupP_Red_AND_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_OXor_SupP_Red_OR(i, r):
        return [f'OXor_SupP_Red_OR_{i}_r{r}_{j}' for j in range(bs)]

    # Initial Degree
    @staticmethod
    def genVars_degree__forward():
        return ['deg_f_' + str(j) for j in range(bs * b)]

    @staticmethod
    def genVars_degree_backward():
        return ['deg_b_' + str(j) for j in range(bs * b)]

    # Match
    @staticmethod
    def genVars_Match_IMC_isWhite(r):
        return [f'Match_IMC_isWhite_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_IP_isWhite(r, pos):
        return [f'Match_IP_isWhite_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_OP_isWhite(r, pos):
        return [f'Match_OP_isWhite_r{r}_{pos}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_OXor_isWhite(r):
        return [f'Match_OXor_isWhite_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Match_Exist(r):
        return [f'Match_Exist_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_Match_Counted(r):
        return [f'G_Match_Counted_r{r}_{j}' for j in range(ColN)]

    # Counter Blue and Red
    @staticmethod
    def genVars_auxi_Xor(i, r):
        return [f'Auxi_Xor_{i}_r{r}_{j}' for j in range(bs)]

    @staticmethod
    def genVars_Blue_Exist(r):
        return [f'Blue_Exist_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_Blue_Counted(r):
        return [f'G_Blue_Counted_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_Red_Exist(r):
        return [f'Red_Exist_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_Red_Counted(r):
        return [f'G_Red_Counted_r{r}_{j}' for j in range(ColN)]


class Constraints_generator:
    def __init__(self, total_round, initial_round, matching_round):
        self.ini_r = initial_round
        self.mat_r = matching_round
        self.TR = total_round
        self.pos = ['L', 'R']
        self.op = ['MC', 'MX', 'Xor']

    def genConstraints_initial_degree(self):
        cons = []
        d1 = Vars_generator.genVars_degree__forward()
        d2 = Vars_generator.genVars_degree_backward()
        IP_1 = Vars_generator.genVars_input_of_round(1, self.ini_r, self.pos[0]) + Vars_generator.genVars_input_of_round(1, self.ini_r, self.pos[1])
        IP_2 = Vars_generator.genVars_input_of_round(2, self.ini_r, self.pos[0]) + Vars_generator.genVars_input_of_round(2, self.ini_r, self.pos[1])
        for bi in range(bs * b):
            cons = cons + [IP_1[bi] + ' + ' + IP_2[bi] + ' >= 1']
            cons = cons + [d1[bi] + ' + ' + IP_2[bi] + ' = 1']
            cons = cons + [d2[bi] + ' + ' + IP_1[bi] + ' = 1']
        return cons

    def genConstraints_round_function(self, r, IMC_1, IMC_2):
        # Operation: MC SB SR MC
        cons = []
        IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[0])
        IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[0])
        IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[0])
        IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[0])
        MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[0])
        Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op[0])
        Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op[0])
        Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op[0])
        # Separate
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                IMC_1[bi],
                IMC_2[bi],
                IMC_SupP_Blue_1[bi],
                IMC_SupP_Blue_2[bi],
                IMC_SupP_Red_1[bi],
                IMC_SupP_Red_2[bi],
                MC_In_isWhite[bi],
                Guess_Blue[bi],
                Guess_Red[bi],
                Guess_Both[bi]
            )
        # MixColumn
        MC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[0])
        MC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[0])
        MC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[0])
        MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[0])
        MC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[0])
        MC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[0])
        MC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[0])
        MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[0])
        OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[0])
        OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[0])
        OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[0])
        OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[0])
        for Coli in range(ColN):
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                column(IMC_SupP_Blue_1, Coli),
                column(IMC_SupP_Blue_2, Coli),
                MC_SupP_Blue_ColExistWhite[Coli],
                MC_SupP_Blue_ColAllGray[Coli],
                column(OMC_SupP_Blue_1, Coli),
                column(OMC_SupP_Blue_2, Coli),
                MC_SupP_Blue_SumGray[Coli],
                MC_CD_Blue[Coli]
            )
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                column(IMC_SupP_Red_1, Coli),
                column(IMC_SupP_Red_2, Coli),
                MC_SupP_Red_ColExistWhite[Coli],
                MC_SupP_Red_ColAllGray[Coli],
                column(OMC_SupP_Red_1, Coli),
                column(OMC_SupP_Red_2, Coli),
                MC_SupP_Red_SumGray[Coli],
                MC_CD_Red[Coli]
            )
        # Merge and ShiftRow
        OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r)
        OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r)
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]], ShiftRow_Inv(OSR_1)[bi])
            cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]], ShiftRow_Inv(OSR_2)[bi])
        # Separate
        IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[1])
        IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[1])
        IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[1])
        IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[1])
        MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[1])
        Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op[1])
        Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op[1])
        Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op[1])
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
        # MixColumn
        MC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[1])
        MC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[1])
        MC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[1])
        MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[1])
        MC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[1])
        MC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[1])
        MC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[1])
        MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[1])
        OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[1])
        OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[1])
        OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[1])
        OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[1])
        for Coli in range(ColN):
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                column(IMC_SupP_Blue_1, Coli),
                column(IMC_SupP_Blue_2, Coli),
                MC_SupP_Blue_ColExistWhite[Coli],
                MC_SupP_Blue_ColAllGray[Coli],
                column(OMC_SupP_Blue_1, Coli),
                column(OMC_SupP_Blue_2, Coli),
                MC_SupP_Blue_SumGray[Coli],
                MC_CD_Blue[Coli]
            )
            cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                column(IMC_SupP_Red_1, Coli),
                column(IMC_SupP_Red_2, Coli),
                MC_SupP_Red_ColExistWhite[Coli],
                MC_SupP_Red_ColAllGray[Coli],
                column(OMC_SupP_Red_1, Coli),
                column(OMC_SupP_Red_2, Coli),
                MC_SupP_Red_SumGray[Coli],
                MC_CD_Red[Coli]
            )
        return cons

    def genConstraints_forward_round(self, r):
        # Operation: SB SR Xor
        cons = []
        IP_L_1 = Vars_generator.genVars_input_of_round(1, r, self.pos[0])
        IP_L_2 = Vars_generator.genVars_input_of_round(2, r, self.pos[0])
        IP_R_1 = Vars_generator.genVars_input_of_round(1, r, self.pos[1])
        IP_R_2 = Vars_generator.genVars_input_of_round(2, r, self.pos[1])
        OP_L_1 = Vars_generator.genVars_output_of_round(1, r, self.pos[0])
        OP_L_2 = Vars_generator.genVars_output_of_round(2, r, self.pos[0])
        OP_R_1 = Vars_generator.genVars_output_of_round(1, r, self.pos[1])
        OP_R_2 = Vars_generator.genVars_output_of_round(2, r, self.pos[1])
        # ShiftRow
        IMC_1 = ShiftRow(IP_L_1)
        IMC_2 = ShiftRow(IP_L_2)
        # Round function
        cons = cons + self.genConstraints_round_function(r, IMC_1, IMC_2)
        # Separate
        IXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_Xor(1, r)
        IXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_Xor(2, r)
        IXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_Xor(1, r)
        IXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_Xor(2, r)
        In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[2])
        Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op[2])
        Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op[2])
        Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op[2])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                IP_R_1[bi],
                IP_R_2[bi],
                IXor_SupP_Blue_1[bi],
                IXor_SupP_Blue_2[bi],
                IXor_SupP_Red_1[bi],
                IXor_SupP_Red_2[bi],
                In_isWhite[bi],
                Guess_Blue[bi],
                Guess_Red[bi],
                Guess_Both[bi]
            )
        # Xor
        OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[1])
        OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[1])
        OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[1])
        OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[1])
        OXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_Xor(1, r)
        OXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_Xor(2, r)
        OXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_Xor(1, r)
        OXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_Xor(2, r)
        CD_Xor_Blue = Vars_generator.genVars_Xor_ConsumedDeg_Blue(r)
        CD_Xor_Red = Vars_generator.genVars_Xor_ConsumedDeg_Red(r)
        OXor_isWhite = Vars_generator.genVars_OXor_isWhite(r)
        OXor_SupP_Blue_AND_1 = Vars_generator.genVars_OXor_SupP_Blue_AND(1, r)
        OXor_SupP_Blue_AND_2 = Vars_generator.genVars_OXor_SupP_Blue_AND(2, r)
        OXor_SupP_Blue_OR_1 = Vars_generator.genVars_OXor_SupP_Blue_OR(1, r)
        OXor_SupP_Blue_OR_2 = Vars_generator.genVars_OXor_SupP_Blue_OR(2, r)
        OXor_SupP_Red_AND_1 = Vars_generator.genVars_OXor_SupP_Red_AND(1, r)
        OXor_SupP_Red_AND_2 = Vars_generator.genVars_OXor_SupP_Red_AND(2, r)
        OXor_SupP_Red_OR_1 = Vars_generator.genVars_OXor_SupP_Red_OR(1, r)
        OXor_SupP_Red_OR_2 = Vars_generator.genVars_OXor_SupP_Red_OR(2, r)
        for bi in range(bs):
            cons = cons + MITMPreConstraints.genConstrains_of_Xor_i(
                IXor_SupP_Blue_1[bi],
                IXor_SupP_Blue_2[bi],
                IXor_SupP_Red_1[bi],
                IXor_SupP_Red_2[bi],
                OMC_SupP_Blue_1[bi],
                OMC_SupP_Blue_2[bi],
                OMC_SupP_Red_1[bi],
                OMC_SupP_Red_2[bi],
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
        # Merge
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[bi], OXor_SupP_Red_1[bi]], OP_L_1[bi])
            cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[bi], OXor_SupP_Red_2[bi]], OP_L_2[bi])
        # Link the output
        cons = cons + MITMPreConstraints.equalConstraints(IMC_1, OP_R_1)
        cons = cons + MITMPreConstraints.equalConstraints(IMC_2, OP_R_2)
        return cons

    def genConstraints_backward_round(self, r):
        # Operation: SB SR Xor
        cons = []
        IP_L_1 = Vars_generator.genVars_input_of_round(1, r, self.pos[0])
        IP_L_2 = Vars_generator.genVars_input_of_round(2, r, self.pos[0])
        IP_R_1 = Vars_generator.genVars_input_of_round(1, r, self.pos[1])
        IP_R_2 = Vars_generator.genVars_input_of_round(2, r, self.pos[1])
        OP_L_1 = Vars_generator.genVars_output_of_round(1, r, self.pos[0])
        OP_L_2 = Vars_generator.genVars_output_of_round(2, r, self.pos[0])
        OP_R_1 = Vars_generator.genVars_output_of_round(1, r, self.pos[1])
        OP_R_2 = Vars_generator.genVars_output_of_round(2, r, self.pos[1])
        # ShiftRow
        cons = cons + MITMPreConstraints.equalConstraints(IP_L_1, ShiftRow_Inv(OP_R_1))
        cons = cons + MITMPreConstraints.equalConstraints(IP_L_2, ShiftRow_Inv(OP_R_2))
        # Round function
        IMC_1 = OP_R_1
        IMC_2 = OP_R_2
        cons = cons + self.genConstraints_round_function(r, IMC_1, IMC_2)
        # Separate
        IXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_Xor(1, r)
        IXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_Xor(2, r)
        IXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_Xor(1, r)
        IXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_Xor(2, r)
        In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[2])
        Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op[2])
        Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op[2])
        Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op[2])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                OP_L_1[bi],
                OP_L_2[bi],
                IXor_SupP_Blue_1[bi],
                IXor_SupP_Blue_2[bi],
                IXor_SupP_Red_1[bi],
                IXor_SupP_Red_2[bi],
                In_isWhite[bi],
                Guess_Blue[bi],
                Guess_Red[bi],
                Guess_Both[bi]
            )
        # Xor
        OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[1])
        OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[1])
        OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[1])
        OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[1])
        OXor_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_Xor(1, r)
        OXor_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_Xor(2, r)
        OXor_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_Xor(1, r)
        OXor_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_Xor(2, r)
        CD_Xor_Blue = Vars_generator.genVars_Xor_ConsumedDeg_Blue(r)
        CD_Xor_Red = Vars_generator.genVars_Xor_ConsumedDeg_Red(r)
        OXor_isWhite = Vars_generator.genVars_OXor_isWhite(r)
        OXor_SupP_Blue_AND_1 = Vars_generator.genVars_OXor_SupP_Blue_AND(1, r)
        OXor_SupP_Blue_AND_2 = Vars_generator.genVars_OXor_SupP_Blue_AND(2, r)
        OXor_SupP_Blue_OR_1 = Vars_generator.genVars_OXor_SupP_Blue_OR(1, r)
        OXor_SupP_Blue_OR_2 = Vars_generator.genVars_OXor_SupP_Blue_OR(2, r)
        OXor_SupP_Red_AND_1 = Vars_generator.genVars_OXor_SupP_Red_AND(1, r)
        OXor_SupP_Red_AND_2 = Vars_generator.genVars_OXor_SupP_Red_AND(2, r)
        OXor_SupP_Red_OR_1 = Vars_generator.genVars_OXor_SupP_Red_OR(1, r)
        OXor_SupP_Red_OR_2 = Vars_generator.genVars_OXor_SupP_Red_OR(2, r)
        for bi in range(bs):
            cons = cons + MITMPreConstraints.genConstrains_of_Xor_i(
                IXor_SupP_Blue_1[bi],
                IXor_SupP_Blue_2[bi],
                IXor_SupP_Red_1[bi],
                IXor_SupP_Red_2[bi],
                OMC_SupP_Blue_1[bi],
                OMC_SupP_Blue_2[bi],
                OMC_SupP_Red_1[bi],
                OMC_SupP_Red_2[bi],
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
        # Merge
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[bi], OXor_SupP_Red_1[bi]], IP_R_1[bi])
            cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[bi], OXor_SupP_Red_2[bi]], IP_R_2[bi])
        return cons

    def genConstraints_link_round(self):
        cons = []
        for r in range(self.TR):
            if r != self.mat_r[0]:
                next_r = (r + 1) % self.TR
                OP_L_1 = Vars_generator.genVars_output_of_round(1, r, self.pos[0])
                OP_L_2 = Vars_generator.genVars_output_of_round(2, r, self.pos[0])
                OP_R_1 = Vars_generator.genVars_output_of_round(1, r, self.pos[1])
                OP_R_2 = Vars_generator.genVars_output_of_round(2, r, self.pos[1])
                IP_L_1 = Vars_generator.genVars_input_of_round(1, next_r, self.pos[0])
                IP_L_2 = Vars_generator.genVars_input_of_round(2, next_r, self.pos[0])
                IP_R_1 = Vars_generator.genVars_input_of_round(1, next_r, self.pos[1])
                IP_R_2 = Vars_generator.genVars_input_of_round(2, next_r, self.pos[1])
                if r == self.TR - 1:
                    cons = cons + MITMPreConstraints.equalConstraints(OP_L_1, IP_R_1)
                    cons = cons + MITMPreConstraints.equalConstraints(OP_L_2, IP_R_2)
                    cons = cons + MITMPreConstraints.equalConstraints(OP_R_1, IP_L_1)
                    cons = cons + MITMPreConstraints.equalConstraints(OP_R_2, IP_L_2)
                else:
                    cons = cons + MITMPreConstraints.equalConstraints(OP_L_1, IP_L_1)
                    cons = cons + MITMPreConstraints.equalConstraints(OP_L_2, IP_L_2)
                    cons = cons + MITMPreConstraints.equalConstraints(OP_R_1, IP_R_1)
                    cons = cons + MITMPreConstraints.equalConstraints(OP_R_2, IP_R_2)
        return cons

    def genConstraints_Match(self):
        cons = []
        IP_L_1 = Vars_generator.genVars_input_of_round(1, self.mat_r[0], self.pos[0])
        IP_L_2 = Vars_generator.genVars_input_of_round(2, self.mat_r[0], self.pos[0])
        IP_R_1 = Vars_generator.genVars_input_of_round(1, self.mat_r[0], self.pos[1])
        IP_R_2 = Vars_generator.genVars_input_of_round(2, self.mat_r[0], self.pos[1])
        OP_L_1 = Vars_generator.genVars_output_of_round(1, self.mat_r[1], self.pos[0])
        OP_L_2 = Vars_generator.genVars_output_of_round(2, self.mat_r[1], self.pos[0])
        OP_R_1 = Vars_generator.genVars_output_of_round(1, self.mat_r[1], self.pos[1])
        OP_R_2 = Vars_generator.genVars_output_of_round(2, self.mat_r[1], self.pos[1])
        for r in self.mat_r:
            IMC_1 = []
            IMC_2 = []
            if r == self.mat_r[0]:
                IMC_1 = ShiftRow(IP_L_1)
                IMC_2 = ShiftRow(IP_L_2)
            if r == self.mat_r[1]:
                IMC_1 = OP_R_1
                IMC_2 = OP_R_2
            # Separate
            IMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(1, r, self.op[0])
            IMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_input_of_MixColumn(2, r, self.op[0])
            IMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(1, r, self.op[0])
            IMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_input_of_MixColumn(2, r, self.op[0])
            MC_In_isWhite = Vars_generator.genVars_In_isWhite(r, self.op[0])
            Guess_Blue = Vars_generator.genVars_Guess_Blue(r, self.op[0])
            Guess_Red = Vars_generator.genVars_Guess_Red(r, self.op[0])
            Guess_Both = Vars_generator.genVars_Guess_Both(r, self.op[0])
            # Separate
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Separate_With_Guess_i(
                    IMC_1[bi],
                    IMC_2[bi],
                    IMC_SupP_Blue_1[bi],
                    IMC_SupP_Blue_2[bi],
                    IMC_SupP_Red_1[bi],
                    IMC_SupP_Red_2[bi],
                    MC_In_isWhite[bi],
                    Guess_Blue[bi],
                    Guess_Red[bi],
                    Guess_Both[bi]
                )
            # MixColumn
            MC_SupP_Blue_ColExistWhite = Vars_generator.genVars_MC_SupP_Blue_ColExistWhite(r, self.op[0])
            MC_SupP_Blue_ColAllGray = Vars_generator.genVars_MC_SupP_Blue_ColAllGray(r, self.op[0])
            MC_SupP_Blue_SumGray = Vars_generator.genVars_MC_SupP_Blue_SumGray(r, self.op[0])
            MC_CD_Blue = Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[0])
            MC_SupP_Red_ColExistWhite = Vars_generator.genVars_MC_SupP_Red_ColExistWhite(r, self.op[0])
            MC_SupP_Red_ColAllGray = Vars_generator.genVars_MC_SupP_Red_ColAllGray(r, self.op[0])
            MC_SupP_Red_SumGray = Vars_generator.genVars_MC_SupP_Red_SumGray(r, self.op[0])
            MC_CD_Red = Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[0])
            OMC_SupP_Blue_1 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(1, r, self.op[0])
            OMC_SupP_Blue_2 = Vars_generator.genVars_SupP_Blue_output_of_MixColumn(2, r, self.op[0])
            OMC_SupP_Red_1 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(1, r, self.op[0])
            OMC_SupP_Red_2 = Vars_generator.genVars_SupP_Red_output_of_MixColumn(2, r, self.op[0])
            for Coli in range(ColN):
                cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
                    column(IMC_SupP_Blue_1, Coli),
                    column(IMC_SupP_Blue_2, Coli),
                    MC_SupP_Blue_ColExistWhite[Coli],
                    MC_SupP_Blue_ColAllGray[Coli],
                    column(OMC_SupP_Blue_1, Coli),
                    column(OMC_SupP_Blue_2, Coli),
                    MC_SupP_Blue_SumGray[Coli],
                    MC_CD_Blue[Coli]
                )
                cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
                    column(IMC_SupP_Red_1, Coli),
                    column(IMC_SupP_Red_2, Coli),
                    MC_SupP_Red_ColExistWhite[Coli],
                    MC_SupP_Red_ColAllGray[Coli],
                    column(OMC_SupP_Red_1, Coli),
                    column(OMC_SupP_Red_2, Coli),
                    MC_SupP_Red_SumGray[Coli],
                    MC_CD_Red[Coli]
                )
            # Merge and ShiftRow
            OSR_1 = Vars_generator.genVars_output_of_ShiftRow(1, r)
            OSR_2 = Vars_generator.genVars_output_of_ShiftRow(2, r)
            for bi in range(bs):
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_1[bi], OMC_SupP_Red_1[bi]], ShiftRow_Inv(OSR_1)[bi])
                cons = cons + MITMPreConstraints.Determine_Allone([OMC_SupP_Blue_2[bi], OMC_SupP_Red_2[bi]], ShiftRow_Inv(OSR_2)[bi])
        # Match
        OSR_matr0_1 = Vars_generator.genVars_output_of_ShiftRow(1, self.mat_r[0])
        OSR_matr0_2 = Vars_generator.genVars_output_of_ShiftRow(2, self.mat_r[0])
        OSR_matr1_1 = Vars_generator.genVars_output_of_ShiftRow(1, self.mat_r[1])
        OSR_matr1_2 = Vars_generator.genVars_output_of_ShiftRow(2, self.mat_r[1])
        IP_L_isWhite = Vars_generator.genVars_Match_IP_isWhite(self.mat_r[0], self.pos[0])
        IP_R_isWhite = Vars_generator.genVars_Match_IP_isWhite(self.mat_r[0], self.pos[1])
        OP_L_isWhite = Vars_generator.genVars_Match_OP_isWhite(self.mat_r[1], self.pos[0])
        OP_R_isWhite = Vars_generator.genVars_Match_OP_isWhite(self.mat_r[1], self.pos[1])
        IMC_matr0_isWhite = Vars_generator.genVars_Match_IMC_isWhite(self.mat_r[0])
        IMC_matr1_isWhite = Vars_generator.genVars_Match_IMC_isWhite(self.mat_r[1])
        OXor_matr0_isWhite = Vars_generator.genVars_Match_OXor_isWhite(self.mat_r[0])
        OXor_matr1_isWhite = Vars_generator.genVars_Match_OXor_isWhite(self.mat_r[1])
        Match_Exist_matr0 = Vars_generator.genVars_Match_Exist(self.mat_r[0])
        Match_Exist_matr1 = Vars_generator.genVars_Match_Exist(self.mat_r[1])
        G_Match_counter_matr0 = Vars_generator.genVars_Match_Counted(self.mat_r[0])
        G_Match_counter_matr1 = Vars_generator.genVars_Match_Counted(self.mat_r[1])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([IP_L_1[bi], IP_L_2[bi]], IP_L_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([IP_R_1[bi], IP_R_2[bi]], IP_R_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OP_L_1[bi], OP_L_2[bi]], OP_L_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OP_R_1[bi], OP_R_2[bi]], OP_R_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OSR_matr0_1[bi], OSR_matr0_2[bi]], IMC_matr0_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_Allzero([OSR_matr1_1[bi], OSR_matr1_2[bi]], IMC_matr1_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_ExistOne([IP_R_isWhite[bi], ShiftRow_Inv(OP_R_isWhite)[bi]], OXor_matr0_isWhite[bi])
        for bi in range(bs):
            cons = cons + MITMPreConstraints.Determine_ExistOne([ShiftRow(IP_L_isWhite)[bi], OP_L_isWhite[bi]], OXor_matr1_isWhite[bi])
        for Coli in range(ColN):
            IMC_matr0_isWhite_Coli = column(IMC_matr0_isWhite, Coli)
            OXor_matr0_isWhite_Coli = column(OXor_matr0_isWhite, Coli)
            cons = cons + [str(RowN + 1) + ' ' + Match_Exist_matr0[Coli] + ' + ' + BasicTools.plusTerms(IMC_matr0_isWhite_Coli + OXor_matr0_isWhite_Coli) + ' <= ' + str(SumIOMC)]
            cons = cons + [str(RowN) + ' ' + Match_Exist_matr0[Coli] + ' + ' + BasicTools.plusTerms(IMC_matr0_isWhite_Coli + OXor_matr0_isWhite_Coli) + ' >= ' + str(RowN)]
            cons = cons + [Match_Exist_matr0[Coli] + ' = 1 -> ' + G_Match_counter_matr0[Coli] + ' + ' + BasicTools.plusTerms(IMC_matr0_isWhite_Coli + OXor_matr0_isWhite_Coli) + ' = ' + str(RowN)]
            cons = cons + [Match_Exist_matr0[Coli] + ' = 0 -> ' + G_Match_counter_matr0[Coli] + ' = 0']
        for Coli in range(ColN):
            IMC_matr1_isWhite_Coli = column(IMC_matr1_isWhite, Coli)
            OXor_matr1_isWhite_Coli = column(OXor_matr1_isWhite, Coli)
            cons = cons + [str(RowN + 1) + ' ' + Match_Exist_matr1[Coli] + ' + ' + BasicTools.plusTerms(IMC_matr1_isWhite_Coli + OXor_matr1_isWhite_Coli) + ' <= ' + str(SumIOMC)]
            cons = cons + [str(RowN) + ' ' + Match_Exist_matr1[Coli] + ' + ' + BasicTools.plusTerms(IMC_matr1_isWhite_Coli + OXor_matr1_isWhite_Coli) + ' >= ' + str(RowN)]
            cons = cons + [Match_Exist_matr1[Coli] + ' = 1 -> ' + G_Match_counter_matr1[Coli] + ' + ' + BasicTools.plusTerms(IMC_matr1_isWhite_Coli + OXor_matr1_isWhite_Coli) + ' = ' + str(RowN)]
            cons = cons + [Match_Exist_matr1[Coli] + ' = 0 -> ' + G_Match_counter_matr1[Coli] + ' = 0']

        # Blue Counter and Red Counter
        # Auxi_Xor_1 = Vars_generator.genVars_auxi_Xor(1, self.mat_r[0])
        # Auxi_Xor_2 = Vars_generator.genVars_auxi_Xor(2, self.mat_r[0])
        # Blue_Exist_matr0 = Vars_generator.genVars_Blue_Exist(self.mat_r[0])
        # Red_Exist_matr0 = Vars_generator.genVars_Red_Exist(self.mat_r[0])
        # Blue_Counter_matr0 = Vars_generator.genVars_Blue_Counted(self.mat_r[0])
        # Red_Counter_matr0 = Vars_generator.genVars_Red_Counted(self.mat_r[0])
        # for bi in range(bs):
        #     cons = cons + MITMPreConstraints.Xor_no_Comsume([IP_R_1[bi], IP_R_2[bi]],
        #                                                     [ShiftRow_Inv(OP_R_1)[bi], ShiftRow_Inv(OP_R_2)[bi]],
        #                                                     [Auxi_Xor_1[bi], Auxi_Xor_2[bi]])
        # for Coli in range(ColN):
        #     IMC_1_Coli = column(OSR_matr0_1, Coli)
        #     OXor_1_Coli = column(Auxi_Xor_1, Coli)
        #     cons = cons + [BasicTools.plusTerms(IMC_1_Coli + OXor_1_Coli) + ' - 5 ' + Blue_Exist_matr0[Coli] + ' >= 0']
        #     cons = cons + [BasicTools.plusTerms(IMC_1_Coli + OXor_1_Coli) + ' - 4 ' + Blue_Exist_matr0[Coli] + ' <= 4']
        #     cons = cons + [Blue_Exist_matr0[Coli] + ' = 1 -> ' + Blue_Counter_matr0[Coli] + ' - ' + BasicTools.minusTerms(IMC_1_Coli + OXor_1_Coli) + ' = -4']
        #     cons = cons + [Blue_Exist_matr0[Coli] + ' = 0 -> ' + Blue_Counter_matr0[Coli] + ' = 0']
        # for Coli in range(ColN):
        #     IMC_2_Coli = column(OSR_matr0_2, Coli)
        #     OXor_2_Coli = column(Auxi_Xor_2, Coli)
        #     cons = cons + [BasicTools.plusTerms(IMC_2_Coli + OXor_2_Coli) + ' - 5 ' + Red_Exist_matr0[Coli] + ' >= 0']
        #     cons = cons + [BasicTools.plusTerms(IMC_2_Coli + OXor_2_Coli) + ' - 4 ' + Red_Exist_matr0[Coli] + ' <= 4']
        #     cons = cons + [Red_Exist_matr0[Coli] + ' = 1 -> ' + Red_Counter_matr0[Coli] + ' - ' + BasicTools.minusTerms(IMC_2_Coli + OXor_2_Coli) + ' = -4']
        #     cons = cons + [Red_Exist_matr0[Coli] + ' = 0 -> ' + Red_Counter_matr0[Coli] + ' = 0']

        # Auxi_Xor_1 = Vars_generator.genVars_auxi_Xor(1, self.mat_r[1])
        # Auxi_Xor_2 = Vars_generator.genVars_auxi_Xor(2, self.mat_r[1])
        # Blue_Exist_matr1 = Vars_generator.genVars_Blue_Exist(self.mat_r[1])
        # Red_Exist_matr1 = Vars_generator.genVars_Red_Exist(self.mat_r[1])
        # Blue_Counter_matr1 = Vars_generator.genVars_Blue_Counted(self.mat_r[1])
        # Red_Counter_matr1 = Vars_generator.genVars_Red_Counted(self.mat_r[1])
        # for bi in range(bs):
        #     cons = cons + MITMPreConstraints.Xor_no_Comsume([ShiftRow(IP_L_1)[bi], ShiftRow(IP_L_2)[bi]],
        #                                                     [OP_L_1[bi], OP_L_2[bi]],
        #                                                     [Auxi_Xor_1[bi], Auxi_Xor_2[bi]])
        # for Coli in range(ColN):
        #     IMC_1_Coli = column(OSR_matr1_1, Coli)
        #     OXor_1_Coli = column(Auxi_Xor_1, Coli)
        #     cons = cons + [BasicTools.plusTerms(IMC_1_Coli + OXor_1_Coli) + ' - 5 ' + Blue_Exist_matr1[Coli] + ' >= 0']
        #     cons = cons + [BasicTools.plusTerms(IMC_1_Coli + OXor_1_Coli) + ' - 4 ' + Blue_Exist_matr1[Coli] + ' <= 4']
        #     cons = cons + [Blue_Exist_matr1[Coli] + ' = 1 -> ' + Blue_Counter_matr1[Coli] + ' - ' + BasicTools.minusTerms(IMC_1_Coli + OXor_1_Coli) + ' = -4']
        #     cons = cons + [Blue_Exist_matr1[Coli] + ' = 0 -> ' + Blue_Counter_matr1[Coli] + ' = 0']
        # for Coli in range(ColN):
        #     IMC_2_Coli = column(OSR_matr1_2, Coli)
        #     OXor_2_Coli = column(Auxi_Xor_2, Coli)
        #     cons = cons + [BasicTools.plusTerms(IMC_2_Coli + OXor_2_Coli) + ' - 5 ' + Red_Exist_matr1[Coli] + ' >= 0']
        #     cons = cons + [BasicTools.plusTerms(IMC_2_Coli + OXor_2_Coli) + ' - 4 ' + Red_Exist_matr1[Coli] + ' <= 4']
        #     cons = cons + [Red_Exist_matr1[Coli] + ' = 1 -> ' + Red_Counter_matr1[Coli] + ' - ' + BasicTools.minusTerms(IMC_2_Coli + OXor_2_Coli) + ' = -4']
        #     cons = cons + [Red_Exist_matr1[Coli] + ' = 0 -> ' + Red_Counter_matr1[Coli] + ' = 0']

        # cons = cons + [BasicTools.plusTerms(Blue_Counter_matr0 + Red_Counter_matr0 + Blue_Counter_matr1 + Red_Counter_matr1) + ' >= 1']
        return cons

    def genConstraints_additional(self):
        cons = []
        CD_Blue = []
        CD_Red = []
        Guess_Blue = []
        Guess_Red = []
        Guess_Both = []
        for r in range(self.TR):
            CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[0])
            CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[0])
            Guess_Blue = Guess_Blue + Vars_generator.genVars_Guess_Blue(r, self.op[0])
            Guess_Red = Guess_Red + Vars_generator.genVars_Guess_Red(r, self.op[0])
            Guess_Both = Guess_Both + Vars_generator.genVars_Guess_Both(r, self.op[0])
            if r not in self.mat_r:
                CD_Blue = CD_Blue + Vars_generator.genVars_MC_ConsumedDeg_Blue(r, self.op[1])
                CD_Blue = CD_Blue + Vars_generator.genVars_Xor_ConsumedDeg_Blue(r)
                CD_Red = CD_Red + Vars_generator.genVars_MC_ConsumedDeg_Red(r, self.op[1])
                CD_Red = CD_Red + Vars_generator.genVars_Xor_ConsumedDeg_Red(r)
                Guess_Blue = Guess_Blue + Vars_generator.genVars_Guess_Blue(r, self.op[1])
                Guess_Blue = Guess_Blue + Vars_generator.genVars_Guess_Blue(r, self.op[2])
                Guess_Red = Guess_Red + Vars_generator.genVars_Guess_Red(r, self.op[1])
                Guess_Red = Guess_Red + Vars_generator.genVars_Guess_Red(r, self.op[2])
                Guess_Both = Guess_Both + Vars_generator.genVars_Guess_Both(r, self.op[1])
                Guess_Both = Guess_Both + Vars_generator.genVars_Guess_Both(r, self.op[2])
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

        G_Match_counter = Vars_generator.genVars_Match_Counted(self.mat_r[0]) + Vars_generator.genVars_Match_Counted(self.mat_r[1])
        GM = 'GMat'
        cons = cons + [GM + ' - ' + BasicTools.minusTerms(G_Match_counter) + ' + ' + BasicTools.plusTerms(Guess_Blue + Guess_Red + Guess_Both) + ' = 0']
        cons = cons + [GM + ' >= 1']
        return cons

    def genConstraints_total(self):
        cons = []
        cons = cons + self.genConstraints_initial_degree()
        if self.ini_r < self.mat_r[0]:
            for r in range(self.ini_r, self.mat_r[0]):
                cons = cons + self.genConstraints_forward_round(r)
            for r in range(0, self.ini_r):
                cons = cons + self.genConstraints_backward_round(r)
            for r in range(self.mat_r[1] + 1, self.TR):
                cons = cons + self.genConstraints_backward_round(r)
        if self.ini_r == self.mat_r[0]:
            for r in range(0, self.ini_r):
                cons = cons + self.genConstraints_backward_round(r)
            for r in range(self.mat_r[1] + 1, self.TR):
                cons = cons + self.genConstraints_backward_round(r)
        if self.ini_r > self.mat_r[1]:
            for r in range(self.ini_r, self.TR):
                cons = cons + self.genConstraints_forward_round(r)
            for r in range(0, self.mat_r[0]):
                cons = cons + self.genConstraints_forward_round(r)
            for r in range(self.mat_r[1] + 1, self.ini_r):
                cons = cons + self.genConstraints_backward_round(r)
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
        cons = cons + ['GObj >= 2']

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
    # TR = 7
    TR = int(sys.argv[1])
    root = f'./Model_Odd/TR{TR}'
    if not os.path.exists(root):
        os.makedirs(root)
    with open(f"./Model_Odd/Result_{TR}.txt", "w") as rd:
        rd.write('TR, ini_r, mat_r: d1, d2, m' + '\n')
        for mat_r in [[i, (i + 1) % TR] for i in range(0, 1)]:
            for ini_r in range(3,4):
                if ini_r != mat_r[1]:
                    filename = f'./Model_Odd/TR{TR}/inir{ini_r}_matr{mat_r}'
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