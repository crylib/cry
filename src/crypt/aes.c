#include <cry/aes.h>
#include "misc.h"
#include <stdint.h>
#include <string.h>


/* Keeps track of key schedule type in the "nr" context element */
#define DECRYPT_STATE       0x80U

/*
 * Forward S-box
 */
static const uint8_t fsbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
    0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
    0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
    0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
    0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
    0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
    0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
    0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
    0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
    0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
    0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
    0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
};

/*
 * Forward tables
 */
#define FTAB \
    V(A5, 63, 63, C6), V(84, 7C, 7C, F8), V(99, 77, 77, EE), V(8D, 7B, 7B, F6), \
    V(0D, F2, F2, FF), V(BD, 6B, 6B, D6), V(B1, 6F, 6F, DE), V(54, C5, C5, 91), \
    V(50, 30, 30, 60), V(03, 01, 01, 02), V(A9, 67, 67, CE), V(7D, 2B, 2B, 56), \
    V(19, FE, FE, E7), V(62, D7, D7, B5), V(E6, AB, AB, 4D), V(9A, 76, 76, EC), \
    V(45, CA, CA, 8F), V(9D, 82, 82, 1F), V(40, C9, C9, 89), V(87, 7D, 7D, FA), \
    V(15, FA, FA, EF), V(EB, 59, 59, B2), V(C9, 47, 47, 8E), V(0B, F0, F0, FB), \
    V(EC, AD, AD, 41), V(67, D4, D4, B3), V(FD, A2, A2, 5F), V(EA, AF, AF, 45), \
    V(BF, 9C, 9C, 23), V(F7, A4, A4, 53), V(96, 72, 72, E4), V(5B, C0, C0, 9B), \
    V(C2, B7, B7, 75), V(1C, FD, FD, E1), V(AE, 93, 93, 3D), V(6A, 26, 26, 4C), \
    V(5A, 36, 36, 6C), V(41, 3F, 3F, 7E), V(02, F7, F7, F5), V(4F, CC, CC, 83), \
    V(5C, 34, 34, 68), V(F4, A5, A5, 51), V(34, E5, E5, D1), V(08, F1, F1, F9), \
    V(93, 71, 71, E2), V(73, D8, D8, AB), V(53, 31, 31, 62), V(3F, 15, 15, 2A), \
    V(0C, 04, 04, 08), V(52, C7, C7, 95), V(65, 23, 23, 46), V(5E, C3, C3, 9D), \
    V(28, 18, 18, 30), V(A1, 96, 96, 37), V(0F, 05, 05, 0A), V(B5, 9A, 9A, 2F), \
    V(09, 07, 07, 0E), V(36, 12, 12, 24), V(9B, 80, 80, 1B), V(3D, E2, E2, DF), \
    V(26, EB, EB, CD), V(69, 27, 27, 4E), V(CD, B2, B2, 7F), V(9F, 75, 75, EA), \
    V(1B, 09, 09, 12), V(9E, 83, 83, 1D), V(74, 2C, 2C, 58), V(2E, 1A, 1A, 34), \
    V(2D, 1B, 1B, 36), V(B2, 6E, 6E, DC), V(EE, 5A, 5A, B4), V(FB, A0, A0, 5B), \
    V(F6, 52, 52, A4), V(4D, 3B, 3B, 76), V(61, D6, D6, B7), V(CE, B3, B3, 7D), \
    V(7B, 29, 29, 52), V(3E, E3, E3, DD), V(71, 2F, 2F, 5E), V(97, 84, 84, 13), \
    V(F5, 53, 53, A6), V(68, D1, D1, B9), V(00, 00, 00, 00), V(2C, ED, ED, C1), \
    V(60, 20, 20, 40), V(1F, FC, FC, E3), V(C8, B1, B1, 79), V(ED, 5B, 5B, B6), \
    V(BE, 6A, 6A, D4), V(46, CB, CB, 8D), V(D9, BE, BE, 67), V(4B, 39, 39, 72), \
    V(DE, 4A, 4A, 94), V(D4, 4C, 4C, 98), V(E8, 58, 58, B0), V(4A, CF, CF, 85), \
    V(6B, D0, D0, BB), V(2A, EF, EF, C5), V(E5, AA, AA, 4F), V(16, FB, FB, ED), \
    V(C5, 43, 43, 86), V(D7, 4D, 4D, 9A), V(55, 33, 33, 66), V(94, 85, 85, 11), \
    V(CF, 45, 45, 8A), V(10, F9, F9, E9), V(06, 02, 02, 04), V(81, 7F, 7F, FE), \
    V(F0, 50, 50, A0), V(44, 3C, 3C, 78), V(BA, 9F, 9F, 25), V(E3, A8, A8, 4B), \
    V(F3, 51, 51, A2), V(FE, A3, A3, 5D), V(C0, 40, 40, 80), V(8A, 8F, 8F, 05), \
    V(AD, 92, 92, 3F), V(BC, 9D, 9D, 21), V(48, 38, 38, 70), V(04, F5, F5, F1), \
    V(DF, BC, BC, 63), V(C1, B6, B6, 77), V(75, DA, DA, AF), V(63, 21, 21, 42), \
    V(30, 10, 10, 20), V(1A, FF, FF, E5), V(0E, F3, F3, FD), V(6D, D2, D2, BF), \
    V(4C, CD, CD, 81), V(14, 0C, 0C, 18), V(35, 13, 13, 26), V(2F, EC, EC, C3), \
    V(E1, 5F, 5F, BE), V(A2, 97, 97, 35), V(CC, 44, 44, 88), V(39, 17, 17, 2E), \
    V(57, C4, C4, 93), V(F2, A7, A7, 55), V(82, 7E, 7E, FC), V(47, 3D, 3D, 7A), \
    V(AC, 64, 64, C8), V(E7, 5D, 5D, BA), V(2B, 19, 19, 32), V(95, 73, 73, E6), \
    V(A0, 60, 60, C0), V(98, 81, 81, 19), V(D1, 4F, 4F, 9E), V(7F, DC, DC, A3), \
    V(66, 22, 22, 44), V(7E, 2A, 2A, 54), V(AB, 90, 90, 3B), V(83, 88, 88, 0B), \
    V(CA, 46, 46, 8C), V(29, EE, EE, C7), V(D3, B8, B8, 6B), V(3C, 14, 14, 28), \
    V(79, DE, DE, A7), V(E2, 5E, 5E, BC), V(1D, 0B, 0B, 16), V(76, DB, DB, AD), \
    V(3B, E0, E0, DB), V(56, 32, 32, 64), V(4E, 3A, 3A, 74), V(1E, 0A, 0A, 14), \
    V(DB, 49, 49, 92), V(0A, 06, 06, 0C), V(6C, 24, 24, 48), V(E4, 5C, 5C, B8), \
    V(5D, C2, C2, 9F), V(6E, D3, D3, BD), V(EF, AC, AC, 43), V(A6, 62, 62, C4), \
    V(A8, 91, 91, 39), V(A4, 95, 95, 31), V(37, E4, E4, D3), V(8B, 79, 79, F2), \
    V(32, E7, E7, D5), V(43, C8, C8, 8B), V(59, 37, 37, 6E), V(B7, 6D, 6D, DA), \
    V(8C, 8D, 8D, 01), V(64, D5, D5, B1), V(D2, 4E, 4E, 9C), V(E0, A9, A9, 49), \
    V(B4, 6C, 6C, D8), V(FA, 56, 56, AC), V(07, F4, F4, F3), V(25, EA, EA, CF), \
    V(AF, 65, 65, CA), V(8E, 7A, 7A, F4), V(E9, AE, AE, 47), V(18, 08, 08, 10), \
    V(D5, BA, BA, 6F), V(88, 78, 78, F0), V(6F, 25, 25, 4A), V(72, 2E, 2E, 5C), \
    V(24, 1C, 1C, 38), V(F1, A6, A6, 57), V(C7, B4, B4, 73), V(51, C6, C6, 97), \
    V(23, E8, E8, CB), V(7C, DD, DD, A1), V(9C, 74, 74, E8), V(21, 1F, 1F, 3E), \
    V(DD, 4B, 4B, 96), V(DC, BD, BD, 61), V(86, 8B, 8B, 0D), V(85, 8A, 8A, 0F), \
    V(90, 70, 70, E0), V(42, 3E, 3E, 7C), V(C4, B5, B5, 71), V(AA, 66, 66, CC), \
    V(D8, 48, 48, 90), V(05, 03, 03, 06), V(01, F6, F6, F7), V(12, 0E, 0E, 1C), \
    V(A3, 61, 61, C2), V(5F, 35, 35, 6A), V(F9, 57, 57, AE), V(D0, B9, B9, 69), \
    V(91, 86, 86, 17), V(58, C1, C1, 99), V(27, 1D, 1D, 3A), V(B9, 9E, 9E, 27), \
    V(38, E1, E1, D9), V(13, F8, F8, EB), V(B3, 98, 98, 2B), V(33, 11, 11, 22), \
    V(BB, 69, 69, D2), V(70, D9, D9, A9), V(89, 8E, 8E, 07), V(A7, 94, 94, 33), \
    V(B6, 9B, 9B, 2D), V(22, 1E, 1E, 3C), V(92, 87, 87, 15), V(20, E9, E9, C9), \
    V(49, CE, CE, 87), V(FF, 55, 55, AA), V(78, 28, 28, 50), V(7A, DF, DF, A5), \
    V(8F, 8C, 8C, 03), V(F8, A1, A1, 59), V(80, 89, 89, 09), V(17, 0D, 0D, 1A), \
    V(DA, BF, BF, 65), V(31, E6, E6, D7), V(C6, 42, 42, 84), V(B8, 68, 68, D0), \
    V(C3, 41, 41, 82), V(B0, 99, 99, 29), V(77, 2D, 2D, 5A), V(11, 0F, 0F, 1E), \
    V(CB, B0, B0, 7B), V(FC, 54, 54, A8), V(D6, BB, BB, 6D), V(3A, 16, 16, 2C)

static const uint32_t ftab[4][256] = {
#define V(a, b, c, d) 0x ## a ## b ## c ## d
    { FTAB },
#undef V
#define V(a, b, c, d) 0x ## b ## c ## d ## a
    { FTAB },
#undef V
#define V(a, b, c, d) 0x ## c ## d ## a ## b
    { FTAB },
#undef V
#define V(a, b, c, d) 0x ## d ## a ## b ## c
    { FTAB }
#undef V
};
#undef FTAB

/*
 * Reverse S-box
 */
static const unsigned char rsbox[256] = {
    0x52, 0x09, 0x6A, 0xD5, 0x30, 0x36, 0xA5, 0x38,
    0xBF, 0x40, 0xA3, 0x9E, 0x81, 0xF3, 0xD7, 0xFB,
    0x7C, 0xE3, 0x39, 0x82, 0x9B, 0x2F, 0xFF, 0x87,
    0x34, 0x8E, 0x43, 0x44, 0xC4, 0xDE, 0xE9, 0xCB,
    0x54, 0x7B, 0x94, 0x32, 0xA6, 0xC2, 0x23, 0x3D,
    0xEE, 0x4C, 0x95, 0x0B, 0x42, 0xFA, 0xC3, 0x4E,
    0x08, 0x2E, 0xA1, 0x66, 0x28, 0xD9, 0x24, 0xB2,
    0x76, 0x5B, 0xA2, 0x49, 0x6D, 0x8B, 0xD1, 0x25,
    0x72, 0xF8, 0xF6, 0x64, 0x86, 0x68, 0x98, 0x16,
    0xD4, 0xA4, 0x5C, 0xCC, 0x5D, 0x65, 0xB6, 0x92,
    0x6C, 0x70, 0x48, 0x50, 0xFD, 0xED, 0xB9, 0xDA,
    0x5E, 0x15, 0x46, 0x57, 0xA7, 0x8D, 0x9D, 0x84,
    0x90, 0xD8, 0xAB, 0x00, 0x8C, 0xBC, 0xD3, 0x0A,
    0xF7, 0xE4, 0x58, 0x05, 0xB8, 0xB3, 0x45, 0x06,
    0xD0, 0x2C, 0x1E, 0x8F, 0xCA, 0x3F, 0x0F, 0x02,
    0xC1, 0xAF, 0xBD, 0x03, 0x01, 0x13, 0x8A, 0x6B,
    0x3A, 0x91, 0x11, 0x41, 0x4F, 0x67, 0xDC, 0xEA,
    0x97, 0xF2, 0xCF, 0xCE, 0xF0, 0xB4, 0xE6, 0x73,
    0x96, 0xAC, 0x74, 0x22, 0xE7, 0xAD, 0x35, 0x85,
    0xE2, 0xF9, 0x37, 0xE8, 0x1C, 0x75, 0xDF, 0x6E,
    0x47, 0xF1, 0x1A, 0x71, 0x1D, 0x29, 0xC5, 0x89,
    0x6F, 0xB7, 0x62, 0x0E, 0xAA, 0x18, 0xBE, 0x1B,
    0xFC, 0x56, 0x3E, 0x4B, 0xC6, 0xD2, 0x79, 0x20,
    0x9A, 0xDB, 0xC0, 0xFE, 0x78, 0xCD, 0x5A, 0xF4,
    0x1F, 0xDD, 0xA8, 0x33, 0x88, 0x07, 0xC7, 0x31,
    0xB1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xEC, 0x5F,
    0x60, 0x51, 0x7F, 0xA9, 0x19, 0xB5, 0x4A, 0x0D,
    0x2D, 0xE5, 0x7A, 0x9F, 0x93, 0xC9, 0x9C, 0xEF,
    0xA0, 0xE0, 0x3B, 0x4D, 0xAE, 0x2A, 0xF5, 0xB0,
    0xC8, 0xEB, 0xBB, 0x3C, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2B, 0x04, 0x7E, 0xBA, 0x77, 0xD6, 0x26,
    0xE1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0C, 0x7D
};

/*
 * Reverse tables
 */
#define RTAB \
    V(50, A7, F4, 51), V(53, 65, 41, 7E), V(C3, A4, 17, 1A), V(96, 5E, 27, 3A), \
    V(CB, 6B, AB, 3B), V(F1, 45, 9D, 1F), V(AB, 58, FA, AC), V(93, 03, E3, 4B), \
    V(55, FA, 30, 20), V(F6, 6D, 76, AD), V(91, 76, CC, 88), V(25, 4C, 02, F5), \
    V(FC, D7, E5, 4F), V(D7, CB, 2A, C5), V(80, 44, 35, 26), V(8F, A3, 62, B5), \
    V(49, 5A, B1, DE), V(67, 1B, BA, 25), V(98, 0E, EA, 45), V(E1, C0, FE, 5D), \
    V(02, 75, 2F, C3), V(12, F0, 4C, 81), V(A3, 97, 46, 8D), V(C6, F9, D3, 6B), \
    V(E7, 5F, 8F, 03), V(95, 9C, 92, 15), V(EB, 7A, 6D, BF), V(DA, 59, 52, 95), \
    V(2D, 83, BE, D4), V(D3, 21, 74, 58), V(29, 69, E0, 49), V(44, C8, C9, 8E), \
    V(6A, 89, C2, 75), V(78, 79, 8E, F4), V(6B, 3E, 58, 99), V(DD, 71, B9, 27), \
    V(B6, 4F, E1, BE), V(17, AD, 88, F0), V(66, AC, 20, C9), V(B4, 3A, CE, 7D), \
    V(18, 4A, DF, 63), V(82, 31, 1A, E5), V(60, 33, 51, 97), V(45, 7F, 53, 62), \
    V(E0, 77, 64, B1), V(84, AE, 6B, BB), V(1C, A0, 81, FE), V(94, 2B, 08, F9), \
    V(58, 68, 48, 70), V(19, FD, 45, 8F), V(87, 6C, DE, 94), V(B7, F8, 7B, 52), \
    V(23, D3, 73, AB), V(E2, 02, 4B, 72), V(57, 8F, 1F, E3), V(2A, AB, 55, 66), \
    V(07, 28, EB, B2), V(03, C2, B5, 2F), V(9A, 7B, C5, 86), V(A5, 08, 37, D3), \
    V(F2, 87, 28, 30), V(B2, A5, BF, 23), V(BA, 6A, 03, 02), V(5C, 82, 16, ED), \
    V(2B, 1C, CF, 8A), V(92, B4, 79, A7), V(F0, F2, 07, F3), V(A1, E2, 69, 4E), \
    V(CD, F4, DA, 65), V(D5, BE, 05, 06), V(1F, 62, 34, D1), V(8A, FE, A6, C4), \
    V(9D, 53, 2E, 34), V(A0, 55, F3, A2), V(32, E1, 8A, 05), V(75, EB, F6, A4), \
    V(39, EC, 83, 0B), V(AA, EF, 60, 40), V(06, 9F, 71, 5E), V(51, 10, 6E, BD), \
    V(F9, 8A, 21, 3E), V(3D, 06, DD, 96), V(AE, 05, 3E, DD), V(46, BD, E6, 4D), \
    V(B5, 8D, 54, 91), V(05, 5D, C4, 71), V(6F, D4, 06, 04), V(FF, 15, 50, 60), \
    V(24, FB, 98, 19), V(97, E9, BD, D6), V(CC, 43, 40, 89), V(77, 9E, D9, 67), \
    V(BD, 42, E8, B0), V(88, 8B, 89, 07), V(38, 5B, 19, E7), V(DB, EE, C8, 79), \
    V(47, 0A, 7C, A1), V(E9, 0F, 42, 7C), V(C9, 1E, 84, F8), V(00, 00, 00, 00), \
    V(83, 86, 80, 09), V(48, ED, 2B, 32), V(AC, 70, 11, 1E), V(4E, 72, 5A, 6C), \
    V(FB, FF, 0E, FD), V(56, 38, 85, 0F), V(1E, D5, AE, 3D), V(27, 39, 2D, 36), \
    V(64, D9, 0F, 0A), V(21, A6, 5C, 68), V(D1, 54, 5B, 9B), V(3A, 2E, 36, 24), \
    V(B1, 67, 0A, 0C), V(0F, E7, 57, 93), V(D2, 96, EE, B4), V(9E, 91, 9B, 1B), \
    V(4F, C5, C0, 80), V(A2, 20, DC, 61), V(69, 4B, 77, 5A), V(16, 1A, 12, 1C), \
    V(0A, BA, 93, E2), V(E5, 2A, A0, C0), V(43, E0, 22, 3C), V(1D, 17, 1B, 12), \
    V(0B, 0D, 09, 0E), V(AD, C7, 8B, F2), V(B9, A8, B6, 2D), V(C8, A9, 1E, 14), \
    V(85, 19, F1, 57), V(4C, 07, 75, AF), V(BB, DD, 99, EE), V(FD, 60, 7F, A3), \
    V(9F, 26, 01, F7), V(BC, F5, 72, 5C), V(C5, 3B, 66, 44), V(34, 7E, FB, 5B), \
    V(76, 29, 43, 8B), V(DC, C6, 23, CB), V(68, FC, ED, B6), V(63, F1, E4, B8), \
    V(CA, DC, 31, D7), V(10, 85, 63, 42), V(40, 22, 97, 13), V(20, 11, C6, 84), \
    V(7D, 24, 4A, 85), V(F8, 3D, BB, D2), V(11, 32, F9, AE), V(6D, A1, 29, C7), \
    V(4B, 2F, 9E, 1D), V(F3, 30, B2, DC), V(EC, 52, 86, 0D), V(D0, E3, C1, 77), \
    V(6C, 16, B3, 2B), V(99, B9, 70, A9), V(FA, 48, 94, 11), V(22, 64, E9, 47), \
    V(C4, 8C, FC, A8), V(1A, 3F, F0, A0), V(D8, 2C, 7D, 56), V(EF, 90, 33, 22), \
    V(C7, 4E, 49, 87), V(C1, D1, 38, D9), V(FE, A2, CA, 8C), V(36, 0B, D4, 98), \
    V(CF, 81, F5, A6), V(28, DE, 7A, A5), V(26, 8E, B7, DA), V(A4, BF, AD, 3F), \
    V(E4, 9D, 3A, 2C), V(0D, 92, 78, 50), V(9B, CC, 5F, 6A), V(62, 46, 7E, 54), \
    V(C2, 13, 8D, F6), V(E8, B8, D8, 90), V(5E, F7, 39, 2E), V(F5, AF, C3, 82), \
    V(BE, 80, 5D, 9F), V(7C, 93, D0, 69), V(A9, 2D, D5, 6F), V(B3, 12, 25, CF), \
    V(3B, 99, AC, C8), V(A7, 7D, 18, 10), V(6E, 63, 9C, E8), V(7B, BB, 3B, DB), \
    V(09, 78, 26, CD), V(F4, 18, 59, 6E), V(01, B7, 9A, EC), V(A8, 9A, 4F, 83), \
    V(65, 6E, 95, E6), V(7E, E6, FF, AA), V(08, CF, BC, 21), V(E6, E8, 15, EF), \
    V(D9, 9B, E7, BA), V(CE, 36, 6F, 4A), V(D4, 09, 9F, EA), V(D6, 7C, B0, 29), \
    V(AF, B2, A4, 31), V(31, 23, 3F, 2A), V(30, 94, A5, C6), V(C0, 66, A2, 35), \
    V(37, BC, 4E, 74), V(A6, CA, 82, FC), V(B0, D0, 90, E0), V(15, D8, A7, 33), \
    V(4A, 98, 04, F1), V(F7, DA, EC, 41), V(0E, 50, CD, 7F), V(2F, F6, 91, 17), \
    V(8D, D6, 4D, 76), V(4D, B0, EF, 43), V(54, 4D, AA, CC), V(DF, 04, 96, E4), \
    V(E3, B5, D1, 9E), V(1B, 88, 6A, 4C), V(B8, 1F, 2C, C1), V(7F, 51, 65, 46), \
    V(04, EA, 5E, 9D), V(5D, 35, 8C, 01), V(73, 74, 87, FA), V(2E, 41, 0B, FB), \
    V(5A, 1D, 67, B3), V(52, D2, DB, 92), V(33, 56, 10, E9), V(13, 47, D6, 6D), \
    V(8C, 61, D7, 9A), V(7A, 0C, A1, 37), V(8E, 14, F8, 59), V(89, 3C, 13, EB), \
    V(EE, 27, A9, CE), V(35, C9, 61, B7), V(ED, E5, 1C, E1), V(3C, B1, 47, 7A), \
    V(59, DF, D2, 9C), V(3F, 73, F2, 55), V(79, CE, 14, 18), V(BF, 37, C7, 73), \
    V(EA, CD, F7, 53), V(5B, AA, FD, 5F), V(14, 6F, 3D, DF), V(86, DB, 44, 78), \
    V(81, F3, AF, CA), V(3E, C4, 68, B9), V(2C, 34, 24, 38), V(5F, 40, A3, C2), \
    V(72, C3, 1D, 16), V(0C, 25, E2, BC), V(8B, 49, 3C, 28), V(41, 95, 0D, FF), \
    V(71, 01, A8, 39), V(DE, B3, 0C, 08), V(9C, E4, B4, D8), V(90, C1, 56, 64), \
    V(61, 84, CB, 7B), V(70, B6, 32, D5), V(74, 5C, 6C, 48), V(42, 57, B8, D0)

static const uint32_t rtab[4][256] = {
#define V(a, b, c, d) 0x ## a ## b ## c ## d
    { RTAB },
#undef V
#define V(a, b, c, d) 0x ## b ## c ## d ## a
    { RTAB },
#undef V
#define V(a, b, c, d) 0x ## c ## d ## a ## b
    { RTAB },
#undef V
#define V(a, b, c, d) 0x ## d ## a ## b ## c
    { RTAB }
#undef V
};
#undef RTAB


static const uint32_t mtab[256] = {
    0x00000000, 0x0b0d090e, 0x161a121c, 0x1d171b12,
    0x2c342438, 0x27392d36, 0x3a2e3624, 0x31233f2a,
    0x58684870, 0x5365417e, 0x4e725a6c, 0x457f5362,
    0x745c6c48, 0x7f516546, 0x62467e54, 0x694b775a,
    0xb0d090e0, 0xbbdd99ee, 0xa6ca82fc, 0xadc78bf2,
    0x9ce4b4d8, 0x97e9bdd6, 0x8afea6c4, 0x81f3afca,
    0xe8b8d890, 0xe3b5d19e, 0xfea2ca8c, 0xf5afc382,
    0xc48cfca8, 0xcf81f5a6, 0xd296eeb4, 0xd99be7ba,
    0x7bbb3bdb, 0x70b632d5, 0x6da129c7, 0x66ac20c9,
    0x578f1fe3, 0x5c8216ed, 0x41950dff, 0x4a9804f1,
    0x23d373ab, 0x28de7aa5, 0x35c961b7, 0x3ec468b9,
    0x0fe75793, 0x04ea5e9d, 0x19fd458f, 0x12f04c81,
    0xcb6bab3b, 0xc066a235, 0xdd71b927, 0xd67cb029,
    0xe75f8f03, 0xec52860d, 0xf1459d1f, 0xfa489411,
    0x9303e34b, 0x980eea45, 0x8519f157, 0x8e14f859,
    0xbf37c773, 0xb43ace7d, 0xa92dd56f, 0xa220dc61,
    0xf66d76ad, 0xfd607fa3, 0xe07764b1, 0xeb7a6dbf,
    0xda595295, 0xd1545b9b, 0xcc434089, 0xc74e4987,
    0xae053edd, 0xa50837d3, 0xb81f2cc1, 0xb31225cf,
    0x82311ae5, 0x893c13eb, 0x942b08f9, 0x9f2601f7,
    0x46bde64d, 0x4db0ef43, 0x50a7f451, 0x5baafd5f,
    0x6a89c275, 0x6184cb7b, 0x7c93d069, 0x779ed967,
    0x1ed5ae3d, 0x15d8a733, 0x08cfbc21, 0x03c2b52f,
    0x32e18a05, 0x39ec830b, 0x24fb9819, 0x2ff69117,
    0x8dd64d76, 0x86db4478, 0x9bcc5f6a, 0x90c15664,
    0xa1e2694e, 0xaaef6040, 0xb7f87b52, 0xbcf5725c,
    0xd5be0506, 0xdeb30c08, 0xc3a4171a, 0xc8a91e14,
    0xf98a213e, 0xf2872830, 0xef903322, 0xe49d3a2c,
    0x3d06dd96, 0x360bd498, 0x2b1ccf8a, 0x2011c684,
    0x1132f9ae, 0x1a3ff0a0, 0x0728ebb2, 0x0c25e2bc,
    0x656e95e6, 0x6e639ce8, 0x737487fa, 0x78798ef4,
    0x495ab1de, 0x4257b8d0, 0x5f40a3c2, 0x544daacc,
    0xf7daec41, 0xfcd7e54f, 0xe1c0fe5d, 0xeacdf753,
    0xdbeec879, 0xd0e3c177, 0xcdf4da65, 0xc6f9d36b,
    0xafb2a431, 0xa4bfad3f, 0xb9a8b62d, 0xb2a5bf23,
    0x83868009, 0x888b8907, 0x959c9215, 0x9e919b1b,
    0x470a7ca1, 0x4c0775af, 0x51106ebd, 0x5a1d67b3,
    0x6b3e5899, 0x60335197, 0x7d244a85, 0x7629438b,
    0x1f6234d1, 0x146f3ddf, 0x097826cd, 0x02752fc3,
    0x335610e9, 0x385b19e7, 0x254c02f5, 0x2e410bfb,
    0x8c61d79a, 0x876cde94, 0x9a7bc586, 0x9176cc88,
    0xa055f3a2, 0xab58faac, 0xb64fe1be, 0xbd42e8b0,
    0xd4099fea, 0xdf0496e4, 0xc2138df6, 0xc91e84f8,
    0xf83dbbd2, 0xf330b2dc, 0xee27a9ce, 0xe52aa0c0,
    0x3cb1477a, 0x37bc4e74, 0x2aab5566, 0x21a65c68,
    0x10856342, 0x1b886a4c, 0x069f715e, 0x0d927850,
    0x64d90f0a, 0x6fd40604, 0x72c31d16, 0x79ce1418,
    0x48ed2b32, 0x43e0223c, 0x5ef7392e, 0x55fa3020,
    0x01b79aec, 0x0aba93e2, 0x17ad88f0, 0x1ca081fe,
    0x2d83bed4, 0x268eb7da, 0x3b99acc8, 0x3094a5c6,
    0x59dfd29c, 0x52d2db92, 0x4fc5c080, 0x44c8c98e,
    0x75ebf6a4, 0x7ee6ffaa, 0x63f1e4b8, 0x68fcedb6,
    0xb1670a0c, 0xba6a0302, 0xa77d1810, 0xac70111e,
    0x9d532e34, 0x965e273a, 0x8b493c28, 0x80443526,
    0xe90f427c, 0xe2024b72, 0xff155060, 0xf418596e,
    0xc53b6644, 0xce366f4a, 0xd3217458, 0xd82c7d56,
    0x7a0ca137, 0x7101a839, 0x6c16b32b, 0x671bba25,
    0x5638850f, 0x5d358c01, 0x40229713, 0x4b2f9e1d,
    0x2264e947, 0x2969e049, 0x347efb5b, 0x3f73f255,
    0x0e50cd7f, 0x055dc471, 0x184adf63, 0x1347d66d,
    0xcadc31d7, 0xc1d138d9, 0xdcc623cb, 0xd7cb2ac5,
    0xe6e815ef, 0xede51ce1, 0xf0f207f3, 0xfbff0efd,
    0x92b479a7, 0x99b970a9, 0x84ae6bbb, 0x8fa362b5,
    0xbe805d9f, 0xb58d5491, 0xa89a4f83, 0xa397468d,
};


#define I0(x) ((x) & 0xff)
#define I1(x) (((x) >> 8) & 0xff)
#define I2(x) (((x) >> 16) & 0xff)
#define I3(x) (((x) >> 24) & 0xff)

#define SBOX(x) \
    ((uint32_t)(fsbox)[I0(x)] | \
    ((uint32_t)(fsbox)[I1(x)] << 8) | \
    ((uint32_t)(fsbox)[I2(x)] << 16) | \
    ((uint32_t)(fsbox)[I3(x)] << 24))

#define ROUND(tab, w0, w1, w2, w3) \
    (tab[0][I0(w0)] ^ \
     tab[1][I1(w1)] ^ \
     tab[2][I2(w2)] ^ \
     tab[3][I3(w3)])

#define FINAL(sbox, w0, w1, w2, w3) \
     ((uint32_t)sbox[I0(w0)] | \
     ((uint32_t)sbox[I1(w1)] << 8) | \
     ((uint32_t)sbox[I2(w2)] << 16) | \
     ((uint32_t)sbox[I3(w3)] << 24))


static void key_invert(cry_aes_ctx *ctx)
{
    unsigned int i, j, k, nr;
    uint32_t key, newkey, tmp;

    nr = ctx->nr & ~DECRYPT_STATE;
    for (i = 0, j = nr * 4; i < j; i += 4, j -= 4) {
        for (k = 0; k < 4; k++)
            CRY_SWAP(ctx->keys[i+k], ctx->keys[j+k]);
    }

    for (i = 4; i < 4 * nr; i++) {
        key = ctx->keys[i];
        newkey = mtab[key & 0xff];
        key >>= 8;
        tmp = mtab[key & 0xff];
        newkey ^= CRY_ROTL32(tmp, 8);
        key >>= 8;
        tmp = mtab[key & 0xff];
        newkey ^= CRY_ROTL32(tmp, 16);
        key >>= 8;
        tmp = mtab[key & 0xff];
        newkey ^= CRY_ROTL32(tmp, 24);
        ctx->keys[i] = newkey;
    }
}

void cry_aes_key_set(cry_aes_ctx *ctx, const unsigned char *key, size_t size)
{
    const unsigned char rcon[10] = {
        0x01, 0x02, 0x04, 0x08, 0x10,
        0x20, 0x40, 0x80, 0x1b, 0x36,
    };
    const unsigned char *rp = rcon;
    unsigned int i, nk, lk;
    uint32_t tmp;

    /* nk = the number of 32 bit words in a key. */
    if (size >= 32)
        nk = 8;
    else if (size >= 24)
        nk = 6;
    else /* size >= 16 */
        nk = 4;
    ctx->nr = nk + 6;
    lk = (CRY_AES_BLOCK_SIZE >> 2) * (ctx->nr + 1);

    /* The first round key is the key itself */
    for (i = 0; i < nk; i++)
        CRY_READ32_LE(ctx->keys[i], key + i*4);

    /* All other round keys are found from the previous round keys */
    for (i = nk; i < lk; i++) {
        tmp = ctx->keys[i-1];
        if (i % nk == 0)
            tmp = SBOX(CRY_ROTL32(tmp, 24)) ^ *rp++;
        else if (nk == 8 && (i%nk) == 4) /* Only for aes-256 */
            tmp = SBOX(tmp);
        ctx->keys[i] = ctx->keys[i-nk] ^ tmp;
    }
}

void cry_aes_encrypt(struct cry_aes_ctx *ctx, unsigned char *dst,
                     const unsigned char *src, size_t size)
{
    uint32_t w0, w1, w2, w3;
    uint32_t t0, t1, t2, t3;
    size_t i;

    if ((ctx->nr & DECRYPT_STATE) != 0) {
        key_invert(ctx);
        ctx->nr &= ~DECRYPT_STATE;
    }

    while (size >= 16) {
        CRY_READ32_LE(w0, src);
        CRY_READ32_LE(w1, src + 4);
        CRY_READ32_LE(w2, src + 8);
        CRY_READ32_LE(w3, src + 12);

        w0 ^= ctx->keys[0];
        w1 ^= ctx->keys[1];
        w2 ^= ctx->keys[2];
        w3 ^= ctx->keys[3];

        for (i = 1; i < ctx->nr; i++) {
            t0 = ROUND(ftab, w0, w1, w2, w3) ^ ctx->keys[4*i];
            t1 = ROUND(ftab, w1, w2, w3, w0) ^ ctx->keys[4*i + 1];
            t2 = ROUND(ftab, w2, w3, w0, w1) ^ ctx->keys[4*i + 2];
            t3 = ROUND(ftab, w3, w0, w1, w2) ^ ctx->keys[4*i + 3];
            w0 = t0;
            w1 = t1;
            w2 = t2;
            w3 = t3;
        }

        t0 = FINAL(fsbox, w0, w1, w2, w3) ^ ctx->keys[4*i];
        t1 = FINAL(fsbox, w1, w2, w3, w0) ^ ctx->keys[4*i + 1];
        t2 = FINAL(fsbox, w2, w3, w0, w1) ^ ctx->keys[4*i + 2];
        t3 = FINAL(fsbox, w3, w0, w1, w2) ^ ctx->keys[4*i + 3];

        CRY_WRITE32_LE(t0, dst);
        CRY_WRITE32_LE(t2, dst + 8);
        CRY_WRITE32_LE(t1, dst + 4);
        CRY_WRITE32_LE(t3, dst + 12);

        size -= CRY_AES_BLOCK_SIZE;
        dst += CRY_AES_BLOCK_SIZE;
        src += CRY_AES_BLOCK_SIZE;
    }
}

void cry_aes_decrypt(struct cry_aes_ctx *ctx, unsigned char *dst,
                     const unsigned char *src, size_t size)
{
    uint32_t w0, w1, w2, w3;
    uint32_t t0, t1, t2, t3;
    size_t i;

    if ((ctx->nr & DECRYPT_STATE) == 0) {
        key_invert(ctx);
        ctx->nr |= DECRYPT_STATE;
    }

    while (size >= 16) {
        CRY_READ32_LE(w0, src);
        CRY_READ32_LE(w1, src + 4);
        CRY_READ32_LE(w2, src + 8);
        CRY_READ32_LE(w3, src + 12);

        w0 ^= ctx->keys[0];
        w1 ^= ctx->keys[1];
        w2 ^= ctx->keys[2];
        w3 ^= ctx->keys[3];

        for (i = 1; i < (ctx->nr & ~DECRYPT_STATE); i++) {
            t0 = ROUND(rtab, w0, w3, w2, w1) ^ ctx->keys[4*i];
            t1 = ROUND(rtab, w1, w0, w3, w2) ^ ctx->keys[4*i + 1];
            t2 = ROUND(rtab, w2, w1, w0, w3) ^ ctx->keys[4*i + 2];
            t3 = ROUND(rtab, w3, w2, w1, w0) ^ ctx->keys[4*i + 3];
            w0 = t0;
            w1 = t1;
            w2 = t2;
            w3 = t3;
        }

        t0 = FINAL(rsbox, w0, w3, w2, w1) ^ ctx->keys[4*i];
        t1 = FINAL(rsbox, w1, w0, w3, w2) ^ ctx->keys[4*i + 1];
        t2 = FINAL(rsbox, w2, w1, w0, w3) ^ ctx->keys[4*i + 2];
        t3 = FINAL(rsbox, w3, w2, w1, w0) ^ ctx->keys[4*i + 3];

        CRY_WRITE32_LE(t0, dst);
        CRY_WRITE32_LE(t2, dst + 8);
        CRY_WRITE32_LE(t1, dst + 4);
        CRY_WRITE32_LE(t3, dst + 12);

        size -= CRY_AES_BLOCK_SIZE;
        dst += CRY_AES_BLOCK_SIZE;
        src += CRY_AES_BLOCK_SIZE;
    }
}

void cry_aes_init(cry_aes_ctx *ctx)
{
    memset(ctx, 0, sizeof(*ctx));
}

void cry_aes_clear(cry_aes_ctx *ctx)
{
    cry_memset(ctx, 0, sizeof(*ctx));
}
