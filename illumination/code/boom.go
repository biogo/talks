// Copyright ©2013 The bíogo.talks Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/biogo/boom"

	"fmt"
	"os"
)

var sampleBam = []byte{
	0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
	0x32, 0x01, 0x7d, 0x50, 0xcb, 0x6e, 0xc2, 0x30, 0x10, 0x0c, 0xc7, 0x20, 0xf1, 0x0b, 0xd5, 0x8a,
	0x4b, 0x4f, 0x76, 0x1e, 0x84, 0x00, 0x39, 0x51, 0x40, 0x6a, 0x2b, 0x41, 0x5a, 0x14, 0xa9, 0x52,
	0x4f, 0x91, 0x63, 0x36, 0x0f, 0x41, 0xec, 0x60, 0x3b, 0xe2, 0x92, 0x8f, 0xaf, 0xd3, 0x0f, 0xe0,
	0xb0, 0xda, 0x99, 0xd1, 0xce, 0x68, 0x35, 0xbb, 0xb7, 0xd3, 0x24, 0x9b, 0x38, 0xce, 0xf6, 0xe3,
	0xe0, 0xfe, 0xa4, 0x49, 0x40, 0x7d, 0x37, 0xfb, 0x4a, 0x7a, 0xa1, 0xa5, 0x32, 0x78, 0x99, 0x4d,
	0xb7, 0xd9, 0xd9, 0xcd, 0xd2, 0xa4, 0x6a, 0x86, 0x20, 0x8e, 0x17, 0x71, 0xb8, 0x8a, 0x82, 0x41,
	0x61, 0x39, 0xa4, 0xfb, 0xdc, 0x0f, 0xfc, 0x70, 0x13, 0xd3, 0x60, 0x70, 0x8f, 0x69, 0xb2, 0x5c,
	0x47, 0xe1, 0x6a, 0xb3, 0xb4, 0x86, 0xef, 0x77, 0xf7, 0xf3, 0x90, 0xec, 0xe4, 0xc3, 0x34, 0x38,
	0x46, 0xfa, 0x34, 0x08, 0xe9, 0xca, 0xdd, 0x1f, 0x93, 0x79, 0xf1, 0x2f, 0x02, 0xd9, 0x03, 0x11,
	0x10, 0x02, 0xb9, 0x41, 0xe4, 0x03, 0x21, 0x05, 0x6a, 0x63, 0x97, 0x51, 0x4d, 0xbb, 0x80, 0x8d,
	0x45, 0x42, 0xb6, 0xec, 0xae, 0x64, 0x2f, 0x2e, 0x96, 0x34, 0xc2, 0x60, 0x85, 0x8a, 0xdc, 0x7b,
	0x76, 0xd3, 0x96, 0x6b, 0xd6, 0x02, 0xe9, 0x20, 0x02, 0x82, 0xb0, 0x1e, 0xed, 0xbc, 0xee, 0xc5,
	0xb5, 0x2d, 0x34, 0x2c, 0x83, 0x10, 0xbc, 0x5a, 0xb6, 0xe8, 0x69, 0x83, 0x5d, 0x8d, 0xa2, 0xf0,
	0x2a, 0xb4, 0x59, 0xa8, 0xbd, 0x53, 0xce, 0x50, 0xf5, 0x55, 0x23, 0xa4, 0x66, 0xde, 0xc9, 0xe2,
	0x9c, 0x03, 0x29, 0x81, 0x7a, 0x5c, 0xde, 0xa4, 0xd2, 0x1d, 0xe3, 0xe8, 0x35, 0x1c, 0x73, 0xae,
	0xf3, 0xf2, 0xd1, 0x8b, 0xbc, 0xe6, 0x94, 0xeb, 0x92, 0x69, 0xc3, 0x80, 0x9c, 0x9f, 0x9e, 0x8d,
	0x6f, 0xc1, 0xa8, 0x8d, 0xa9, 0x91, 0x9f, 0xb7, 0xac, 0xa3, 0xf6, 0xc5, 0xf9, 0x6c, 0x6a, 0x5b,
	0x75, 0x5e, 0xec, 0x3c, 0xed, 0xce, 0xb9, 0xbe, 0xfe, 0x3a, 0x7f, 0xe6, 0x6e, 0x80, 0x2d, 0x85,
	0x01, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42,
	0x43, 0x02, 0x00, 0xf0, 0x01, 0x8d, 0xd3, 0xcb, 0x8e, 0xda, 0x30, 0x14, 0x06, 0x60, 0x43, 0x77,
	0xa3, 0x19, 0x71, 0x55, 0xc9, 0xc0, 0x90, 0x1b, 0x81, 0x24, 0x04, 0x27, 0x71, 0x12, 0x3b, 0xb6,
	0x21, 0x80, 0x03, 0xaa, 0xb2, 0x61, 0xd6, 0x55, 0x55, 0x29, 0xab, 0x76, 0xd1, 0x6e, 0xfa, 0x08,
	0x05, 0xb1, 0x98, 0x87, 0xea, 0x13, 0xf4, 0xa5, 0xa6, 0x1e, 0x15, 0x54, 0x21, 0x44, 0xd4, 0x23,
	0x59, 0xb2, 0xe4, 0x85, 0x3f, 0xfd, 0xe7, 0x9c, 0x3d, 0xf8, 0x5b, 0x9f, 0x7f, 0x3f, 0x81, 0xe6,
	0xeb, 0xaf, 0x76, 0x4d, 0xde, 0x6d, 0x79, 0x5e, 0x4f, 0xf5, 0xf6, 0x86, 0x70, 0x89, 0xa2, 0x28,
	0x2e, 0x11, 0x8e, 0x71, 0xf9, 0x21, 0x06, 0x3f, 0xea, 0x00, 0x1c, 0x0f, 0x85, 0x91, 0x1f, 0x8a,
	0xc2, 0x3a, 0x76, 0x95, 0xe2, 0xb8, 0x7d, 0x69, 0x15, 0xa6, 0xce, 0x3d, 0x38, 0x5f, 0xcf, 0x54,
	0x9f, 0x2d, 0x89, 0xa6, 0xc0, 0xcc, 0x53, 0xfd, 0x2c, 0x83, 0xc6, 0xd8, 0x0e, 0xb9, 0xd5, 0xb0,
	0xa2, 0xa9, 0x61, 0x4d, 0xa6, 0xbd, 0x87, 0x80, 0x7c, 0x14, 0x9b, 0xfa, 0x6e, 0xfb, 0x29, 0x66,
	0xe0, 0x79, 0xb7, 0x01, 0x9b, 0xdd, 0xa6, 0xfe, 0xf3, 0x84, 0xf8, 0x22, 0x11, 0x8d, 0x1b, 0x08,
	0x56, 0xa2, 0x18, 0xbf, 0x19, 0x50, 0x78, 0x36, 0x38, 0x96, 0x54, 0x08, 0x71, 0x7c, 0xd9, 0x16,
	0xfb, 0xee, 0xc1, 0xcc, 0x5b, 0x56, 0x57, 0x17, 0x29, 0xcd, 0xf3, 0x0c, 0xf3, 0x3c, 0xf7, 0x46,
	0x09, 0x5f, 0x05, 0x5e, 0x96, 0x2d, 0x22, 0x9b, 0x88, 0x35, 0x4d, 0x16, 0x22, 0x77, 0xc6, 0xab,
	0xf5, 0xca, 0xd2, 0x33, 0x21, 0x0d, 0xe0, 0xc2, 0x00, 0xce, 0x86, 0xaf, 0x15, 0x06, 0x5a, 0xa2,
	0x14, 0x4b, 0x48, 0x82, 0xa3, 0xb3, 0xe1, 0x32, 0x02, 0xcb, 0x68, 0x8b, 0x7c, 0x8d, 0xb1, 0x58,
	0x8b, 0xd5, 0x32, 0xcb, 0xb2, 0x94, 0xf0, 0x39, 0x83, 0x66, 0xef, 0xce, 0x62, 0xc1, 0xf8, 0xd1,
	0x8e, 0xbd, 0x78, 0xd4, 0x34, 0x71, 0xf4, 0x00, 0x3a, 0xae, 0x01, 0x6e, 0xe7, 0x50, 0x65, 0x40,
	0x32, 0x03, 0x4a, 0xd3, 0x32, 0x49, 0xc2, 0x2a, 0x03, 0x41, 0x9a, 0x86, 0xed, 0xa1, 0x01, 0x99,
	0x3b, 0x54, 0x8c, 0xd0, 0x68, 0xde, 0x77, 0x6c, 0xc7, 0x35, 0x07, 0xfe, 0xa4, 0x5e, 0xb7, 0x91,
	0xd6, 0x6e, 0x0d, 0xa6, 0xe3, 0x51, 0x5f, 0x91, 0x86, 0xda, 0x85, 0xa1, 0xf6, 0x5f, 0x06, 0x39,
	0x10, 0x54, 0xfe, 0x8f, 0x13, 0x54, 0x65, 0xd0, 0x9d, 0xd9, 0x84, 0x21, 0x88, 0xd9, 0x9c, 0x60,
	0x4f, 0x31, 0xa7, 0xae, 0x63, 0xf8, 0x94, 0x07, 0x7d, 0x75, 0xd0, 0x19, 0xc7, 0x89, 0x3b, 0x34,
	0xbd, 0xfe, 0xbd, 0x1d, 0xda, 0xb7, 0x7b, 0xf1, 0xad, 0xda, 0x90, 0xa4, 0x91, 0x74, 0xc8, 0x40,
	0x4e, 0x86, 0x2b, 0x41, 0xde, 0xd2, 0x67, 0x90, 0xb2, 0x65, 0x82, 0x28, 0x31, 0xdf, 0x35, 0x83,
	0x49, 0xaf, 0xad, 0x52, 0xbb, 0xd1, 0x18, 0x42, 0x0c, 0xef, 0xde, 0x53, 0xd8, 0x7d, 0x34, 0x9f,
	0xd4, 0x81, 0xa2, 0x7a, 0xda, 0x75, 0x2f, 0xf6, 0x27, 0xc3, 0xf7, 0xaa, 0xc5, 0x40, 0x72, 0x31,
	0x52, 0xd9, 0x10, 0xc6, 0xc8, 0xbf, 0x20, 0x2e, 0x27, 0x52, 0x42, 0xf4, 0xd8, 0x5f, 0x08, 0x81,
	0x63, 0xc6, 0x97, 0x2c, 0x74, 0xc3, 0x2c, 0x98, 0xcd, 0x39, 0x27, 0x96, 0x8e, 0xb8, 0xe7, 0x10,
	0xce, 0xd3, 0x34, 0x22, 0x0b, 0xa8, 0xf5, 0xfa, 0xd7, 0x41, 0xfc, 0x01, 0x00, 0xa5, 0x23, 0x08,
	0x9e, 0x03, 0x00, 0x00, 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
	0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
}

func init() {
	f, err := os.Create("sample.bam")
	if err != nil {
		panic(err)
	}
	defer f.Close()
	_, err = f.Write(sampleBam)
	if err != nil {
		panic(err)
	}
}

func main() {
	bf, err := boom.OpenBAM("sample.bam")
	if err != nil {
		fmt.Fprintf(os.Stderr, "could not open file: %v\n", err)
		os.Exit(1)
	}
	fmt.Println(bf.RefNames())

	for {
		r, _, err := bf.Read()
		if err != nil {
			break
		}
		fmt.Println(r)
	}
}
