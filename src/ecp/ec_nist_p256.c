/*
 * Copyright (c) 2013-2016, Davide Galassi. All rights reserved.
 *
 * This file is part of CRY software.
 *
 * CRY is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with CRY; if not, see <http://www.gnu/licenses/>.
 */

#include "cry/ec.h"

/* NIST Curve P-256 parameters */

int cry_ec_set_nist_p256(cry_ec *ec)
{
    int res;

    unsigned char p[] = {
        0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x01,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF,
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
    };

    unsigned char b[] = {
        0x5A, 0xC6, 0x35, 0xD8, 0xAA, 0x3A, 0x93, 0xE7,
        0xB3, 0xEB, 0xBD, 0x55, 0x76, 0x98, 0x86, 0xBC,
        0x65, 0x1D, 0x06, 0xB0, 0xCC, 0x53, 0xB0, 0xF6,
        0x3B, 0xCE, 0x3C, 0x3E, 0x27, 0xD2, 0x60, 0x4B
    };

    unsigned char n[] = {
        0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x00, 0x00, 0x00,
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
        0xBC, 0xE6, 0xFA, 0xAD, 0xA7, 0x17, 0x9E, 0x84,
        0xF3, 0xB9, 0xCA, 0xC2, 0xFC, 0x63, 0x25, 0x51
    };

    unsigned char gx[] = {
        0x6B, 0x17, 0xD1, 0xF2, 0xE1, 0x2C, 0x42, 0x47,
        0xF8, 0xBC, 0xE6, 0xE5, 0x63, 0xA4, 0x40, 0xF2,
        0x77, 0x03, 0x7D, 0x81, 0x2D, 0xEB, 0x33, 0xA0,
        0xF4, 0xA1, 0x39, 0x45, 0xD8, 0x98, 0xC2, 0x96
    };

    unsigned char gy[] = {
        0x4F, 0xE3, 0x42, 0xE2, 0xFE, 0x1A, 0x7F, 0x9B,
        0x8E, 0xE7, 0xEB, 0x4A, 0x7C, 0x0F, 0x9E, 0x16,
        0x2B, 0xCE, 0x33, 0x57, 0x6B, 0x31, 0x5E, 0xCE,
        0xCB, 0xB6, 0x40, 0x68, 0x37, 0xBF, 0x51, 0xF5
    };

    if ((res = cry_mpi_load_bin(&ec->p, p, sizeof(p))) != 0)
        return res;
    if ((res = cry_mpi_set_int(&ec->a, -3)) != 0)
        return res;
    if ((res = cry_mpi_load_bin(&ec->b, b, sizeof(b))) != 0)
        return res;
    if ((res = cry_mpi_load_bin(&ec->g.x, gx, sizeof(gx))) != 0)
        return res;
    if ((res = cry_mpi_load_bin(&ec->g.y, gy, sizeof(gy))) != 0)
        return res;
    if ((res = cry_mpi_load_bin(&ec->n, n, sizeof(n))) != 0)
        return res;
    return 0;
}

int cry_ec_init_nist_p256(cry_ec *ec)
{
    int res;

    if ((res = cry_ec_init(ec)) == 0) {
        if ((res = cry_ec_set_nist_p256(ec)) != 0)
            cry_ec_clear(ec);
    }
    return res;
}

