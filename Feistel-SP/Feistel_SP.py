

from gurobipy import *
from utils import *
import sys
import os


class Vars_generate:
    @staticmethod
    def genVars_output_of_MixColumn(r, i, op=''):
        return [f'OMC_{op}_{i}_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_SupP_Blue_output_of_MixColumn(r, i, op=''):
        return [f'OMC_SupP_Blue_{op}_{i}_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_SupP_Red_output_of_MixColumn(r, i, op=''):
        return [f'OMC_SupP_Red_{op}_{i}_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_output_of_SubBytes(r, i):
        return [f'OSB_{i}_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_SupP_Blue_output_of_SubBytes(r, i):
        return [f'OSB_SupP_Blue_{i}_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_SupP_Red_output_of_SubBytes(r, i):
        return [f'OSB_SupP_Red_{i}_r{r}_{j}' for j in range(RowN)]

    # For Separate
    @staticmethod
    def genVars_IMC_isWhite(r):
        return [f'IMC_isWhite_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_IXor_isWhite(r):
        return [f'IXor_isWhite_r{r}_{j}' for j in range(RowN)]

    # For MixColumn
    @staticmethod
    def genVars_MC_SupP_Blue_ColExistWhite(r):
        return [f'MC_SupP_Blue_ColExistWhite_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Blue_ColAllGray(r):
        return [f'MC_SupP_Blue_ColAllGray_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Blue_SumGray(r):
        return [f'G_SupP_Blue_SumGray_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_ConsumedDeg_Blue(r):
        return [f'G_CD_MC_Blue_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Red_ColExistWhite(r):
        return [f'MC_SupP_Red_ColExistWhite_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Red_ColAllGray(r):
        return [f'MC_SupP_Red_ColAllGray_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_SupP_Red_SumGray(r):
        return [f'G_SupP_Red_SumGray_r{r}_{j}' for j in range(ColN)]

    @staticmethod
    def genVars_MC_ConsumedDeg_Red(r):
        return [f'G_CD_MC_Red_r{r}_{j}' for j in range(ColN)]

    # For Xor
    @staticmethod
    def genVars_OXor_SupP_Blue_AND(r, i):
        return [f'OXor_SupP_Blue_AND_r{r}_{i}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_OXor_SupP_Blue_OR(r, i):
        return [f'OXor_SupP_Blue_OR_r{r}_{i}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_OXor_SupP_Red_AND(r, i):
        return [f'OXor_SupP_Red_AND_r{r}_{i}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_OXor_SupP_Red_OR(r, i):
        return [f'OXor_SupP_Red_OR_r{r}_{i}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_Xor_ConsumedDeg_Blue(r):
        return [f'CD_Xor_Blue_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_Xor_ConsumedDeg_Red(r):
        return [f'CD_Xor_Red_r{r}_{j}' for j in range(RowN)]

    @staticmethod
    def genVars_OXor_isWhite(r):
        return [f'OXor_isWhite_r{r}_{j}' for j in range(RowN)]

    # Initial Degree
    @staticmethod
    def genVars_ini_degree_Forward():
        return [f'deg_f_{j}' for j in range(2 * RowN)]

    @staticmethod
    def genVars_ini_degree_Backward():
        return [f'deg_b_{j}' for j in range(2 * RowN)]

    # Match
    @staticmethod
    def genVars_AllZero(r):
        return [f'AZ_{r}_{i}' for i in range(RowN)]

    @staticmethod
    def genVars_ExistDoubleZero():
        return [[f'EDZ_{i}_{j}' for j in range(RowN)] for i in range(4)]

    @staticmethod
    def genVars_Match():
        return [[f'm_{i}_{j}' for j in range(RowN)] for i in range(2)]


class Constraints_generate:
    def __init__(self, TR, ini_r, mat_r):
        self.TR = TR
        self.ini_r = ini_r
        self.mat_r = mat_r

    def genConstraints_ini_degree(self):
        cons = []
        d1 = Vars_generate.genVars_ini_degree_Forward()
        d2 = Vars_generate.genVars_ini_degree_Backward()
        OMC_U_1 = Vars_generate.genVars_output_of_MixColumn(self.ini_r[0], 1)
        OMC_U_2 = Vars_generate.genVars_output_of_MixColumn(self.ini_r[0], 2)
        OMC_L_1 = Vars_generate.genVars_output_of_MixColumn(self.ini_r[1], 1)
        OMC_L_2 = Vars_generate.genVars_output_of_MixColumn(self.ini_r[1], 2)

        for i in range(RowN):
            cons = cons + [OMC_U_1[i] + ' + ' + OMC_U_2[i] + ' >= 1']
            cons = cons + [d1[i] + ' + ' + OMC_U_2[i] + ' = 1']
            cons = cons + [d2[i] + ' + ' + OMC_U_1[i] + ' = 1']
        for i in range(RowN):
            cons = cons + [OMC_L_1[i] + ' + ' + OMC_L_2[i] + ' >= 1']
            cons = cons + [d1[RowN+ i] + ' + ' + OMC_L_2[i] + ' = 1']
            cons = cons + [d2[RowN + i] + ' + ' + OMC_L_1[i] + ' = 1']

        return cons

    @staticmethod
    def genConstrains_of_MixColumn_Xor(input1_r, input2_r, output_r):
        cons = []

        IXor_1 = Vars_generate.genVars_output_of_MixColumn(input1_r, 1)
        IXor_2 = Vars_generate.genVars_output_of_MixColumn(input1_r, 2)
        IMC_1 = Vars_generate.genVars_output_of_SubBytes(input2_r, 1)
        IMC_2 = Vars_generate.genVars_output_of_SubBytes(input2_r, 2)
        OMC_next_r_1 = Vars_generate.genVars_output_of_MixColumn(output_r, 1)
        OMC_next_r_2 = Vars_generate.genVars_output_of_MixColumn(output_r, 2)

        # Separate for Input
        IXor_SupP_Blue_1 = Vars_generate.genVars_SupP_Blue_output_of_MixColumn(input1_r, 1)
        IXor_SupP_Blue_2 = Vars_generate.genVars_SupP_Blue_output_of_MixColumn(input1_r, 2)
        IXor_SupP_Red_1 = Vars_generate.genVars_SupP_Red_output_of_MixColumn(input1_r, 1)
        IXor_SupP_Red_2 = Vars_generate.genVars_SupP_Red_output_of_MixColumn(input1_r, 2)
        IXor_In_isWhite = Vars_generate.genVars_IXor_isWhite(input1_r)
        IMC_SupP_Blue_1 = Vars_generate.genVars_SupP_Blue_output_of_SubBytes(input2_r, 1)
        IMC_SupP_Blue_2 = Vars_generate.genVars_SupP_Blue_output_of_SubBytes(input2_r, 2)
        IMC_SupP_Red_1 = Vars_generate.genVars_SupP_Red_output_of_SubBytes(input2_r, 1)
        IMC_SupP_Red_2 = Vars_generate.genVars_SupP_Red_output_of_SubBytes(input2_r, 2)
        IMC_In_isWhite = Vars_generate.genVars_IMC_isWhite(input2_r)

        # Variables for MixColumn
        op = 'MX'
        OMC_SupP_Blue_1 = Vars_generate.genVars_SupP_Blue_output_of_MixColumn(output_r, 1, op)
        OMC_SupP_Blue_2 = Vars_generate.genVars_SupP_Blue_output_of_MixColumn(output_r, 2, op)
        MC_SupP_Blue_ColExistWhite = Vars_generate.genVars_MC_SupP_Blue_ColExistWhite(output_r)[0]
        MC_SupP_Blue_ColAllGray = Vars_generate.genVars_MC_SupP_Blue_ColAllGray(output_r)[0]
        G_SupP_Blue_SumGray = Vars_generate.genVars_MC_SupP_Blue_SumGray(output_r)[0]
        G_CD_MC_Blue = Vars_generate.genVars_MC_ConsumedDeg_Blue(output_r)[0]
        OMC_SupP_Red_1 = Vars_generate.genVars_SupP_Red_output_of_MixColumn(output_r, 1, op)
        OMC_SupP_Red_2 = Vars_generate.genVars_SupP_Red_output_of_MixColumn(output_r, 2, op)
        MC_SupP_Red_ColExistWhite = Vars_generate.genVars_MC_SupP_Red_ColExistWhite(output_r)[0]
        MC_SupP_Red_ColAllGray = Vars_generate.genVars_MC_SupP_Red_ColAllGray(output_r)[0]
        G_SupP_Red_SumGray = Vars_generate.genVars_MC_SupP_Red_SumGray(output_r)[0]
        G_CD_MC_Red = Vars_generate.genVars_MC_ConsumedDeg_Red(output_r)[0]

        # Output of Xor
        OXor_SupP_Blue_1 = Vars_generate.genVars_SupP_Blue_output_of_MixColumn(output_r, 1)
        OXor_SupP_Blue_2 = Vars_generate.genVars_SupP_Blue_output_of_MixColumn(output_r, 2)
        OXor_SupP_Red_1 = Vars_generate.genVars_SupP_Red_output_of_MixColumn(output_r, 1)
        OXor_SupP_Red_2 = Vars_generate.genVars_SupP_Red_output_of_MixColumn(output_r, 2)
        OXor_SupP_Blue_AND_1 = Vars_generate.genVars_OXor_SupP_Blue_AND(output_r, 1)
        OXor_SupP_Blue_AND_2 = Vars_generate.genVars_OXor_SupP_Blue_AND(output_r, 2)
        OXor_SupP_Blue_OR_1 = Vars_generate.genVars_OXor_SupP_Blue_OR(output_r, 1)
        OXor_SupP_Blue_OR_2 = Vars_generate.genVars_OXor_SupP_Blue_OR(output_r, 2)
        OXor_SupP_Red_AND_1 = Vars_generate.genVars_OXor_SupP_Red_AND(output_r, 1)
        OXor_SupP_Red_AND_2 = Vars_generate.genVars_OXor_SupP_Red_AND(output_r, 2)
        OXor_SupP_Red_OR_1 = Vars_generate.genVars_OXor_SupP_Red_OR(output_r, 1)
        OXor_SupP_Red_OR_2 = Vars_generate.genVars_OXor_SupP_Red_OR(output_r, 2)
        G_CD_Xor_Blue = Vars_generate.genVars_Xor_ConsumedDeg_Blue(output_r)
        G_CD_Xor_Red = Vars_generate.genVars_Xor_ConsumedDeg_Red(output_r)
        OXor_isWhite = Vars_generate.genVars_OXor_isWhite(output_r)

        for Rowi in range(RowN):
            cons = cons + MITMPreConstraints.Separate_without_Guess_i(
                IXor_1[Rowi],
                IXor_2[Rowi],
                IXor_SupP_Blue_1[Rowi],
                IXor_SupP_Blue_2[Rowi],
                IXor_SupP_Red_1[Rowi],
                IXor_SupP_Red_2[Rowi],
                IXor_In_isWhite[Rowi]
            )
            cons = cons + MITMPreConstraints.Separate_without_Guess_i(
                IMC_1[Rowi],
                IMC_2[Rowi],
                IMC_SupP_Blue_1[Rowi],
                IMC_SupP_Blue_2[Rowi],
                IMC_SupP_Red_1[Rowi],
                IMC_SupP_Red_2[Rowi],
                IMC_In_isWhite[Rowi]
            )
        cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Blue(
            IMC_SupP_Blue_1,
            IMC_SupP_Blue_2,
            MC_SupP_Blue_ColExistWhite,
            MC_SupP_Blue_ColAllGray,
            OMC_SupP_Blue_1,
            OMC_SupP_Blue_2,
            G_SupP_Blue_SumGray,
            G_CD_MC_Blue
        )
        cons = cons + MITMPreConstraints.genSubConstraints_MC_SupP__Red(
            IMC_SupP_Red_1,
            IMC_SupP_Red_2,
            MC_SupP_Red_ColExistWhite,
            MC_SupP_Red_ColAllGray,
            OMC_SupP_Red_1,
            OMC_SupP_Red_2,
            G_SupP_Red_SumGray,
            G_CD_MC_Red
        )

        for Rowi in range(RowN):
            cons = cons + MITMPreConstraints.genConstrains_of_Xor_i(
                IXor_SupP_Blue_1[Rowi],
                IXor_SupP_Blue_2[Rowi],
                IXor_SupP_Red_1[Rowi],
                IXor_SupP_Red_2[Rowi],
                OMC_SupP_Blue_1[Rowi],
                OMC_SupP_Blue_2[Rowi],
                OMC_SupP_Red_1[Rowi],
                OMC_SupP_Red_2[Rowi],
                OXor_SupP_Blue_1[Rowi],
                OXor_SupP_Blue_2[Rowi],
                OXor_SupP_Red_1[Rowi],
                OXor_SupP_Red_2[Rowi],
                G_CD_Xor_Blue[Rowi],
                G_CD_Xor_Red[Rowi],
                OXor_isWhite[Rowi],
                OXor_SupP_Blue_AND_1[Rowi],
                OXor_SupP_Blue_AND_2[Rowi],
                OXor_SupP_Blue_OR_1[Rowi],
                OXor_SupP_Blue_OR_2[Rowi],
                OXor_SupP_Red_AND_1[Rowi],
                OXor_SupP_Red_AND_2[Rowi],
                OXor_SupP_Red_OR_1[Rowi],
                OXor_SupP_Red_OR_2[Rowi]
            )

        for Rowi in range(RowN):
            cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_1[Rowi], OXor_SupP_Red_1[Rowi]], OMC_next_r_1[Rowi])
            cons = cons + MITMPreConstraints.Determine_Allone([OXor_SupP_Blue_2[Rowi], OXor_SupP_Red_2[Rowi]], OMC_next_r_2[Rowi])

        return cons

    @staticmethod
    def genConstrains_of_equal(r):
        cons = []
        OMC1 = Vars_generate.genVars_output_of_MixColumn(r, 1)
        OMC2 = Vars_generate.genVars_output_of_MixColumn(r, 2)
        OSB1 = Vars_generate.genVars_output_of_SubBytes(r, 1)
        OSB2 = Vars_generate.genVars_output_of_SubBytes(r, 2)
        cons = cons + MITMPreConstraints.equalConstraints(OMC1, OSB1)
        cons = cons + MITMPreConstraints.equalConstraints(OMC2, OSB2)
        return cons

    def genConstrains_of_last_round(self, n):
        cons = []
        for i in range(n):
            OMC_1 = Vars_generate.genVars_output_of_MixColumn(i, 1)
            OMC_2 = Vars_generate.genVars_output_of_MixColumn(i, 2)
            OMC_last_r_1 = Vars_generate.genVars_output_of_MixColumn(self.TR - 1 - i, 1)
            OMC_last_r_2 = Vars_generate.genVars_output_of_MixColumn(self.TR - 1 - i, 2)
            cons = cons + MITMPreConstraints.equalConstraints(OMC_1, OMC_last_r_1)
            cons = cons + MITMPreConstraints.equalConstraints(OMC_2, OMC_last_r_2)
        return cons

    def genConstraints_Match(self):
        cons = []
        M = Vars_generate.genVars_Match()
        EDZ = Vars_generate.genVars_ExistDoubleZero()
        AZ_FOR_EDZ = [[[] for j in range(RowN)] for i in range(4)]
        R_FOR_EDZ = [[] for i in range(4)]

        for r in range(self.TR):
            OSB1 = Vars_generate.genVars_output_of_SubBytes(r, 1)
            OSB2 = Vars_generate.genVars_output_of_SubBytes(r, 2)
            AllZero = Vars_generate.genVars_AllZero(r)
            for i in range(RowN):
                cons = cons + MITMPreConstraints.Determine_Allzero([OSB1[i], OSB2[i]], AllZero[i])

        if self.mat_r >= self.ini_r[1]:
            # Forward
            if self.mat_r == self.TR - 1:
                R_FOR_EDZ[0] = [r for r in range(self.ini_r[1], self.mat_r + 1) if (self.mat_r - r) % 2 == 1]
                R_FOR_EDZ[1] = [r for r in range(self.ini_r[1], self.mat_r + 1) if (self.mat_r - r) % 2 == 0]
            else:
                R_FOR_EDZ[0] = [r for r in range(self.ini_r[1], self.mat_r + 1) if (self.mat_r - r) % 2 == 0]
                R_FOR_EDZ[1] = [r for r in range(self.ini_r[1], self.mat_r + 1) if (self.mat_r - r) % 2 == 1]
            # Backward
            if self.mat_r == self.TR - 1:
                R_FOR_EDZ[2] = [r for r in range(0, self.ini_r[1]) if r % 2 == 1]
                R_FOR_EDZ[3] = [r for r in range(0, self.ini_r[1]) if r % 2 == 0]
            elif (self.TR - 1 - self.mat_r) % 2 == 0:
                R_FOR_EDZ[2] = [r for r in range(0, self.ini_r[1]) if r % 2 == 0]
                R_FOR_EDZ[3] = [r for r in range(0, self.ini_r[1]) if r % 2 == 1]
            else:
                R_FOR_EDZ[2] = [r for r in range(0, self.ini_r[1]) if r % 2 == 1]
                R_FOR_EDZ[3] = [r for r in range(0, self.ini_r[1]) if r % 2 == 0]
            R_FOR_EDZ[2] = R_FOR_EDZ[2] + [r for r in range(self.mat_r + 1, self.TR) if (r - self.mat_r) % 2 == 0]
            R_FOR_EDZ[3] = R_FOR_EDZ[3] + [r for r in range(self.mat_r + 1, self.TR) if (r - self.mat_r) % 2 == 1]

        elif self.mat_r < self.ini_r[0]:
            # Forward
            if self.mat_r % 2 == 0:
                R_FOR_EDZ[0] = [r for r in range(self.ini_r[1], self.TR) if (self.TR - 1 - r) % 2 == 0]
                R_FOR_EDZ[1] = [r for r in range(self.ini_r[1], self.TR) if (self.TR - 1 - r) % 2 == 1]
            else:
                R_FOR_EDZ[0] = [r for r in range(self.ini_r[1], self.TR) if (self.TR - 1 - r) % 2 == 1]
                R_FOR_EDZ[1] = [r for r in range(self.ini_r[1], self.TR) if (self.TR - 1 - r) % 2 == 0]

            R_FOR_EDZ[0] = R_FOR_EDZ[0] + [r for r in range(0, self.mat_r + 1) if (self.mat_r - r) % 2 == 0]
            R_FOR_EDZ[1] = R_FOR_EDZ[1] + [r for r in range(0, self.mat_r + 1) if (self.mat_r - r) % 2 == 1]
            # Backward
            R_FOR_EDZ[2] = [r for r in range(self.mat_r + 1, self.ini_r[1]) if (r - self.mat_r) % 2 == 0]
            R_FOR_EDZ[3] = [r for r in range(self.mat_r + 1, self.ini_r[1]) if (r - self.mat_r) % 2 == 1]

        for i in range(4):
            for j in range(RowN):
                for r in R_FOR_EDZ[i]:
                    AZ_FOR_EDZ[i][j].append(Vars_generate.genVars_AllZero(r)[j])
        for i in range(4):
            for j in range(RowN):
                cons = cons + BasicTools.N_OR_(AZ_FOR_EDZ[i][j], EDZ[i][j])

        GMat = 'GMat'
        for i in range(b):
            for j in range(RowN):
                cons = cons + BasicTools.AND([EDZ[i][j], EDZ[i + 2][j]], M[i][j])
        cons = cons + [GMat + ' - ' + BasicTools.minusTerms(M[0] + M[1]) + ' = 0']
        cons = cons + ['GMat >= 1']
        # return R_FOR_EDZ
        return cons

    def genConstraints_of_additional(self):
        cons = []
        CD_Blue = []
        CD_Red = []
        if self.mat_r < self.ini_r[0]:
            for r in range(0, self.TR + 1):
                if r not in self.ini_r:
                    CD_Blue = CD_Blue + Vars_generate.genVars_MC_ConsumedDeg_Blue(r) + Vars_generate.genVars_Xor_ConsumedDeg_Blue(r)
                    CD_Red = CD_Red  + Vars_generate.genVars_MC_ConsumedDeg_Red(r) + Vars_generate.genVars_Xor_ConsumedDeg_Red(r)
        if self.mat_r >= self.ini_r[1]:
            for r in range(-1, self.TR):
                if r not in self.ini_r:
                    CD_Blue = CD_Blue + Vars_generate.genVars_MC_ConsumedDeg_Blue(r) + Vars_generate.genVars_Xor_ConsumedDeg_Blue(r)
                    CD_Red = CD_Red  + Vars_generate.genVars_MC_ConsumedDeg_Red(r) + Vars_generate.genVars_Xor_ConsumedDeg_Red(r)

        d1 = Vars_generate.genVars_ini_degree_Forward()
        d2 = Vars_generate.genVars_ini_degree_Backward()

        Deg1 = 'GDeg1'
        Deg2 = 'GDeg2'

        if len(CD_Blue) > 0:
            cons = cons + [Deg1 + ' - ' + BasicTools.minusTerms(d1) + ' + ' + BasicTools._plusTerms(CD_Blue) + ' = 0']
        else:
            cons = cons + [Deg1 + ' - ' + BasicTools.minusTerms(d1) + ' = 0']

        if len(CD_Red) > 0:
            cons = cons + [Deg2 + ' - ' + BasicTools.minusTerms(d2) + ' + ' + BasicTools._plusTerms(CD_Red) + ' = 0']
        else:
            cons = cons + [Deg2 + ' - ' + BasicTools.minusTerms(d2) + ' = 0']

        cons = cons + [Deg1 + ' >= 1']
        cons = cons + [Deg2 + ' >= 1']

        return cons

    def genConstraints_of_total(self):
        cons = []
        cons = cons + self.genConstraints_ini_degree()

        for r in range(self.TR):
            cons = cons + self.genConstrains_of_equal(r)

        if self.mat_r < self.ini_r[0]:
            # Forward
            for r in range(self.ini_r[1] + 1, self.TR):
                cons = cons + self.genConstrains_of_MixColumn_Xor(r - 2, r - 1, r)
            # Forward
            if self.mat_r == 0:
                cons = cons + self.genConstrains_of_last_round(1)
            elif self.mat_r == 1:
                cons = cons + self.genConstrains_of_last_round(2)
            elif self.mat_r >= 2:
                cons = cons + self.genConstrains_of_last_round(3)
                for r in range(3, self.mat_r + 1):
                    cons = cons + self.genConstrains_of_MixColumn_Xor(r - 2, r - 1, r)
            # Backward
            for r in range(self.mat_r + 1, self.ini_r[0]):
                cons = cons + self.genConstrains_of_MixColumn_Xor(r + 2, r + 1, r)
        if self.mat_r >= self.ini_r[1]:
            # Forward
            for r in range(self.ini_r[1] + 1, self.mat_r + 1):
                cons = cons + self.genConstrains_of_MixColumn_Xor(r - 2, r - 1, r)
            # Backward
            for r in range(0, self.ini_r[0]):
                cons = cons + self.genConstrains_of_MixColumn_Xor(r + 2, r + 1, r)
            # Backward
            if self.mat_r == self.TR - 2:
                cons = cons + self.genConstrains_of_last_round(1)
            elif self.mat_r == self.TR - 3:
                cons = cons + self.genConstrains_of_last_round(2)
            elif self.mat_r <= self.TR - 4:
                cons = cons + self.genConstrains_of_last_round(3)
                for r in range(self.mat_r + 1, self.TR - 3):
                    cons = cons + self.genConstrains_of_MixColumn_Xor(r + 2, r + 1, r)

        cons = cons + self.genConstraints_Match()
        cons = cons + self.genConstraints_of_additional()

        return cons

    def genModel(self, filename):
        V = set([])
        cons = []
        cons = cons + self.genConstraints_of_total()

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
    TR = int(sys.argv[1])
    Dir = f'./Model/TR{TR}'
    Result_File = f"./Model/Result_{TR}.txt"
    if not os.path.exists(Dir):
        os.mkdir(Dir)
    with open(Result_File, "w") as rd:
        rd.write('TR, ini_r, mat_r: d1, d2, m' + '\n')
        for ini_r in [[i, i + 1] for i in range(1, TR - 2)]:
            for mat_r in range(0, TR):
                if mat_r != ini_r[0]:
                    filename = f'./Model/TR{TR}/inir{ini_r}_matr{mat_r}'
                    A = Constraints_generate(TR, ini_r, mat_r)
                    A.genModel(filename)
                    Model = read(filename + '.lp')
                    # Model.setParam('TimeLimit', 120 * 60)
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
                                temp = temp.replace('-', ' ')
                                temp = temp.split()
                                Sol[temp[0]] = int(temp[1])
                        rd.write(str(TR) + ',' + str(ini_r) + ',' + str(mat_r) + ':')
                        rd.write(str(Sol['GDeg1']) + ',' + str(Sol['GDeg2']) + ',' + str(Sol['GMat']) + '\n')
                        rd.flush()




