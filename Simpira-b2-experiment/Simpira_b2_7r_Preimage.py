
import time
import copy
import random

bs = 16
b = 2

sbox = [
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b,
    0xfe, 0xd7, 0xab, 0x76, 0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, 0xb7, 0xfd, 0x93, 0x26,
    0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2,
    0xeb, 0x27, 0xb2, 0x75, 0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, 0x53, 0xd1, 0x00, 0xed,
    0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f,
    0x50, 0x3c, 0x9f, 0xa8, 0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, 0xcd, 0x0c, 0x13, 0xec,
    0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14,
    0xde, 0x5e, 0x0b, 0xdb, 0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, 0xe7, 0xc8, 0x37, 0x6d,
    0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f,
    0x4b, 0xbd, 0x8b, 0x8a, 0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, 0xe1, 0xf8, 0x98, 0x11,
    0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f,
    0xb0, 0x54, 0xbb, 0x16
]

inv_sbox = [
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e,
    0x81, 0xf3, 0xd7, 0xfb, 0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87,
    0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb, 0x54, 0x7b, 0x94, 0x32,
    0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49,
    0x6d, 0x8b, 0xd1, 0x25, 0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16,
    0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92, 0x6c, 0x70, 0x48, 0x50,
    0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05,
    0xb8, 0xb3, 0x45, 0x06, 0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02,
    0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b, 0x3a, 0x91, 0x11, 0x41,
    0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8,
    0x1c, 0x75, 0xdf, 0x6e, 0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89,
    0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b, 0xfc, 0x56, 0x3e, 0x4b,
    0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59,
    0x27, 0x80, 0xec, 0x5f, 0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d,
    0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef, 0xa0, 0xe0, 0x3b, 0x4d,
    0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63,
    0x55, 0x21, 0x0c, 0x7d
]


# ShiftRow and SubBytes
def sub_bytes(state):
    res = [0] * bs
    shift = [0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11]
    for i in range(bs):
        res[i] = sbox[state[shift[i]]]
    return res


# multiplication in the finite field
def ff_mult(a, b):
    p = 0
    tmp = 0
    for i in range(8):
        if b & 1 == 1:
            p ^= a
        tmp = a & 0x80
        a <<= 1
        if tmp == 0x80:
            a ^= 0x1b
        b >>= 1
    return p % 256


# MixColumns
def mix_columns(state):
    res = [0] * bs
    for i in range(4):
        res[0 + 4 * i] = (ff_mult(state[0 + 4 * i], 2)
                          ^ ff_mult(state[3 + 4 * i], 1)
                          ^ ff_mult(state[2 + 4 * i], 1)
                          ^ ff_mult(state[1 + 4 * i], 3))
        res[1 + 4 * i] = (ff_mult(state[1 + 4 * i], 2)
                          ^ ff_mult(state[0 + 4 * i], 1)
                          ^ ff_mult(state[3 + 4 * i], 1)
                          ^ ff_mult(state[2 + 4 * i], 3))
        res[2 + 4 * i] = (ff_mult(state[2 + 4 * i], 2)
                          ^ ff_mult(state[1 + 4 * i], 1)
                          ^ ff_mult(state[0 + 4 * i], 1)
                          ^ ff_mult(state[3 + 4 * i], 3))
        res[3 + 4 * i] = (ff_mult(state[3 + 4 * i], 2)
                          ^ ff_mult(state[2 + 4 * i], 1)
                          ^ ff_mult(state[1 + 4 * i], 1)
                          ^ ff_mult(state[0 + 4 * i], 3))
    return res


# Inv MixColumns
def inv_mix_columns(state):
    res = [0] * bs
    for i in range(4):
        res[0 + 4 * i] = (ff_mult(state[0 + 4 * i], 0xe)
                          ^ ff_mult(state[3 + 4 * i], 0x9)
                          ^ ff_mult(state[2 + 4 * i], 0xd)
                          ^ ff_mult(state[1 + 4 * i], 0xb))
        res[1 + 4 * i] = (ff_mult(state[1 + 4 * i], 0xe)
                          ^ ff_mult(state[0 + 4 * i], 0x9)
                          ^ ff_mult(state[3 + 4 * i], 0xd)
                          ^ ff_mult(state[2 + 4 * i], 0xb))
        res[2 + 4 * i] = (ff_mult(state[2 + 4 * i], 0xe)
                          ^ ff_mult(state[1 + 4 * i], 0x9)
                          ^ ff_mult(state[0 + 4 * i], 0xd)
                          ^ ff_mult(state[3 + 4 * i], 0xb))
        res[3 + 4 * i] = (ff_mult(state[3 + 4 * i], 0xe)
                          ^ ff_mult(state[2 + 4 * i], 0x9)
                          ^ ff_mult(state[1 + 4 * i], 0xd)
                          ^ ff_mult(state[0 + 4 * i], 0xb))
    return res


# aes round
def aes_round(state, round_constant):
    res = [0] * bs
    # subbytes and shiftrows
    shift = [0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11]
    for i in range(bs):
        res[i] = sbox[state[shift[i]]]
    # mix columns
    res = mix_columns(res)
    for i in range(len(res)):
        res[i] = res[i] ^ round_constant[i]
    return res


# Round constant depending on b (number of branches) and i (current index)
def const(i, b=2):
    tmp0 = ((i ^ b) % 2**8)
    tmp1 = ((i ^ b) // 2**8) % (2**8)
    tmp2 = ((i ^ b) // 2**16) % (2**8)
    tmp3 = ((i ^ b) // 2**24) % (2**8)
    return [
        0x00 ^ tmp0, tmp1, tmp2, tmp3, 0x10 ^ tmp0, tmp1, tmp2, tmp3,
        0x20 ^ tmp0, tmp1, tmp2, tmp3, 0x30 ^ tmp0, tmp1, tmp2, tmp3
    ]


# XOR of two states
def xor(x, y):
    return [x[i] ^ y[i] for i in range(len(x))]


# Round function pi
def pi(x, i, b=2):
    return aes_round(aes_round(x, const(i, b=b)), [0] * 16)


# Inverse of pi
def invpi(x, i, b=2):
    return inv_aes_round(inv_aes_round(x, [0] * 16), const(i, b=b))


# Full or reduced Simpira-2
def simpira2(x, nrounds=15):
    for r in range(nrounds):
        x[(r + 1) % 2] = xor(x[(r + 1) % 2], pi(x[r % 2], r+1, b=2))


def simpira2_MitM_7r_start_at_round4(x):
    # - Forward round 4, 5, 6
    x_forward = copy.deepcopy(x)
    for r in range(4, 7):
        x_forward[(r + 1) % 2] = xor(x_forward[(r + 1) % 2], pi(x_forward[r % 2], r+1, b=2))
    # - Backward round 3, 2, 1, 0
    x_backward = copy.deepcopy(x)
    for r in range(3, -1, -1):
        x_backward[(r + 1) % 2] = xor(x_backward[(r + 1) % 2], pi(x_backward[r % 2], (r + 1), b=2))
    target = [[0 for bi in range(bs)], [0 for bi in range(bs)]]
    for bi in range(bs):
        target[0][bi] = x_forward[0][bi] ^ x_backward[0][bi]
        target[1][bi] = x_forward[1][bi] ^ x_backward[1][bi]
    target_inv_MixColumn = [inv_mix_columns(target[0]), inv_mix_columns(target[1])]
    if target_inv_MixColumn[0][9] == 0 and target_inv_MixColumn[0][12] == 0:
        print('input_inv:  ', aes_state_to_hex(inv_mix_columns(x_backward[0])), aes_state_to_hex(inv_mix_columns(x_backward[1])))
        print('output_inv: ', aes_state_to_hex(inv_mix_columns(x_forward[0])), aes_state_to_hex(inv_mix_columns(x_forward[1])))
        print('target_inv: ', aes_state_to_hex(target_inv_MixColumn[0]))
        simpira2(x_backward, 7)
        if x_backward == x_forward:
            return 1
    else:

        return 0


# conversion of an AES state into 4 32-bit hex numbers representing the columns
def aes_state_to_hex(state):
    # map state to columns
    res = [0, 0, 0, 0]
    for i in range(4):
        for j in range(4):
            res[i] += state[4 * i + j] * ((256)**(3 - j))
    return [hex(t) for t in res]


def test():
    # - input_inv_mc:  90d64cee 5dceafc3 c0600c7b 1a4ecd95 cbce2e53 fe452225 e49464ea 31d57501
    # - output_inv_mc: ce623383 274f3cb0 bf603c92 1a43ae10 ca060030 b89b5a75 4352d9a3 fb5c6f95
    input_state_inv_MC = [
        [0x90, 0xd6, 0x4c, 0xee, 0x5d, 0xce, 0xaf, 0xc3, 0xc0, 0x60, 0x0c, 0x7b, 0x1a, 0x4e, 0xcd, 0x95],
        [0xcb, 0xce, 0x2e, 0x53, 0xfe, 0x45, 0x22, 0x25, 0xe4, 0x94, 0x64, 0xea, 0x31, 0xd5, 0x75, 0x01]
    ]
    input_state = [mix_columns(input_state_inv_MC[0]), mix_columns(input_state_inv_MC[1])]
    output_state = [[0] * 16, [0] * 16]
    for i in range(b):
        for bi in range(bs):
            output_state[i][bi] = input_state[i][bi]
    simpira2(output_state, 7)
    output_state_inv_MC = [inv_mix_columns(output_state[0]), inv_mix_columns(output_state[1])]
    assert aes_state_to_hex(output_state_inv_MC[0])[0] == '0xce623383'

    target_inv_MC = [[0] * 16, [0] * 16]
    for i in range(b):
        for bi in range(bs):
            target_inv_MC[i][bi] = input_state_inv_MC[i][bi] ^ output_state_inv_MC[i][bi]
    print(aes_state_to_hex(target_inv_MC[0])[2])
    print(aes_state_to_hex(target_inv_MC[0])[3])


def compute_partial_preimage():
    red_position = [0, 7, 10, 13]
    blue_gray_position = [1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 14, 15]
    counter = 0
    start = time.time()
    for i in range(2 ** 10):
        print('episode: ', i)
        start_state1 = [random.randint(0, 255) for bi in range(bs)]  # A^3_MC1
        start_state1[5] = 0
        start_state1[14] = 0
        g_r = [random.randint(0, 255) for bi in range(bs)]
        L = {}
        for red_val in range(256):
            # - Construct the starting states
            start_state1[5] = red_val # B^4
            A_3_SR2 = sub_bytes(aes_round(start_state1, const(4)))
            # - Fix the blue and gray cells to be 0
            A_3_SR2_copy = copy.deepcopy(A_3_SR2)
            for pos in blue_gray_position:
                A_3_SR2_copy[pos] = 0
            start_state2 = xor(mix_columns(A_3_SR2_copy), g_r) # A^4
            # - Forward
            A_5 = xor(start_state1, pi(start_state2, 5))
            A_5_SR2 = sub_bytes(aes_round(A_5, const(6)))
            tmp = hex(((A_3_SR2[9] ^ A_5_SR2[9]) * 256) ^ A_3_SR2[12] ^ A_5_SR2[12])
            if tmp in L:
                L[tmp].append([start_state2, red_val])
            else:
                L[tmp] = []
                L[tmp].append([start_state2, red_val])
        
        for blue_val in range(256):
            start_state1[14] = blue_val  # B^4
            A_3_SR2 = sub_bytes(aes_round(start_state1, const(4)))
            for pos in red_position:
                A_3_SR2[pos] = 0
            B_3 = xor(mix_columns(A_3_SR2), g_r)
            A_3 = start_state1
            # - Backward
            A_2 = B_3
            B_2 = xor(A_3, pi(A_2, 3))
            A_1_SR2 = sub_bytes(aes_round(B_2, const(2)))
            tmp = hex((A_1_SR2[9] * 256) ^ A_1_SR2[12])
            if tmp in L:
                for v in L[tmp]:
                    start_state1[5] = v[1]
                    if simpira2_MitM_7r_start_at_round4([v[0], start_state1]):
                        counter += 1
    end = time.time()
    print(end - start)
    print(counter)


if __name__ == '__main__':
    compute_partial_preimage()
    # test()
