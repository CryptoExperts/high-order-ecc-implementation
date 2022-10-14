// -------------------------------------------------------------------------- //
//---- Language:         C
//---- Module name:      ecdh.c
//---- Developed:        Sonia Belaïd & Matthieu Rivain
//---- Contact:          {sonia.belaid,matthieu.rivain}@cryptoexperts.com
//---- Copyright:        CryptoExperts
// -------------------------------------------------------------------------- //

#include "lowlevel.h"
#include "ecdh.h"

// -------------------------------------------------------------------------- //
// ---                       GLOBAL VARIABLES                             --- //
// -------------------------------------------------------------------------- //

/** random buffer in case of external source of randomness */
unsigned int ECDH_RANDOM_BUFFER[128];
/** counter for random integers in the random buffer */
unsigned int delta;
/** working memory zone */
unsigned int workzone[1000];
/** number of needed registers depending on the countermeasure */
unsigned char nb_reg;
/** pointers for low lovel operations */
unsigned int *in1, *in2, aux, *out; 

/******************************************************************************/
/***                               MACRO                                    ***/
/******************************************************************************/

/** activation of the countermeasure */
#define ECC_CM
/** type */
#define u64 unsigned long long int
/** sizes */
#define OPSIZE(x)	x
#define MODSIZE(x)	(x+2)
#define LEN 		(getL())
/** number of registers, depending on the countermeasure */
#define NBACC		((unsigned char) 1 << (nb_reg - 2))
/** memory workzone */
#define MOD			workzone
#define AFF_X		(MOD + MODSIZE(LEN) + 3*OPSIZE(LEN))
#define AFF_Y		(AFF_X + OPSIZE(LEN))
#define JAC_X		(AFF_Y + 3*OPSIZE(LEN))
#define JAC_Y		(JAC_X + OPSIZE(LEN))
#define JAC_Z		(JAC_Y + OPSIZE(LEN))
#define S_MASK		(JAC_Z + OPSIZE(LEN))	
#define SCAL		(S_MASK + 12*(OPSIZE(LEN)+1))					
#define SCAL1		(SCAL + OPSIZE(LEN)+1)				
#define SCAL2		(SCAL1 + OPSIZE(LEN)+1)				
#define SCAL3		(SCAL2 + OPSIZE(LEN)+1)				
#define EXPX		(SCAL3 + OPSIZE(LEN)+1)				
#define EXPY		(EXPX + OPSIZE(LEN))				
#define EXPZ		(EXPY + OPSIZE(LEN))				
#define EXPW		(EXPZ + OPSIZE(LEN))				
#define ACCX		(EXPW + OPSIZE(LEN))				
#define ACCY		(ACCX + NBACC*(OPSIZE(LEN)))		
#define ACCZ		(ACCY + NBACC*(OPSIZE(LEN)))		
#define TMP1		(ACCZ + NBACC*(OPSIZE(LEN)))		
#define TMP2		(TMP1 + OPSIZE(LEN))				
#define TMP3		(TMP2 + OPSIZE(LEN))				
#define TMP4		(TMP3 + OPSIZE(LEN))				
#define TMP5		(TMP4 + OPSIZE(LEN))				
#define TMP6		(TMP5 + OPSIZE(LEN))				
#define TMP7		(TMP6 + OPSIZE(LEN))
#define TMPX		SCAL
#define TMPY		(SCAL + OPSIZE(LEN))
#define TMPZ		(TMPY+ OPSIZE(LEN))
#define RESX		(SCAL + OPSIZE(LEN))
#define RESY		(RESX+ OPSIZE(LEN))

/******************************************************************************/
/******************************************************************************/
/***                         ENTRY POINTS                                   ***/
/******************************************************************************/
/******************************************************************************/

// -------------------------------------------------------------------------- //
// ---                      Set Random Buffer                             --- //
// -------------------------------------------------------------------------- //

void set_random_buffer(unsigned int* random_bytes, unsigned int len)
{
	int i;
	
	for (i=0; i<len; i++)
	{
		ECDH_RANDOM_BUFFER[i] = random_bytes[i];
	}
}


/******************************************************************************/
/******************************************************************************/
/***                  ECDH						                            ***/
/******************************************************************************/
/******************************************************************************/

// -------------------------------------------------------------------------- //
// ---                    ROUTINE DECLARATIONS                            --- //
// -------------------------------------------------------------------------- //

/**
 * ECDH primitive (scalar multiplication)
 * - level:
 * 		if level=1: addition of a multiple (32-bit random) of the order to the scalar
 * 		otherwise (level=2): Euclidean division of the scalar with 64 bits
 * 							 hence 3 scalar multiplication
 * - output: output array (result of the scalar multiplication)
 * - modulus: modulus
 * - data_a: domain parameter
 * - data_b: domain parameter
 * - data_r: curve order
 * - cm_param: masking order
 * - scalar: scalar
 * - x_coord: point X-coordinate
 * - y_coord: point Y-coordinate
 * - len: modulus length
 */
void ecdh(int level, unsigned int* output, unsigned int* modulus, unsigned int* data_a, unsigned int* data_b, unsigned int* data_r, unsigned int* cm_param, unsigned int* scalar, unsigned int* x_coord, unsigned int* y_coord, int len);
/**
 * First steps of the scalar multiplication:
 * 		Modulus setup and Montgomery precomputations
 * - modulus
 * - len: modulus length
 */
void ecdh_scalar_multiplication_init(unsigned int* modulus, int len);
/**
 * Core of the scalar multiplication:
 * - level:
 * 		if level=1: addition of a multiple (32-bit random) of the order to the scalar
 * 		otherwise (level=2): Euclidean division of the scalar with 64 bits
 * 							 hence 3 scalar multiplication
 * - param: masking order
 * - data_a: domain parameter
 * - data_r: curve order
 * - scalar: scalar
 * - x_coord: point X-coordinate
 * - y_coord: point Y-coordinate
 */
void ecdh_scalar_multiplication_core(int level, unsigned char param, unsigned int* data_a, unsigned int* data_r, unsigned int* scalar, unsigned int* x_coord, unsigned int* y_coord);
/**
 * Scalar multiplication with COZ coordinates
 * - param: masking order
 * - data_a: domain parameter
 * - data_r: curve order
 * - scalar: scalar
 * - x_coord: point X-coordinate
 * - y_coord: point Y-coordinate
 */
void ecdh_scalar_point_mult_binary_coz(unsigned char param, unsigned int* data_a, unsigned int* data_r, unsigned int *scalar, unsigned int scalar_len, unsigned int* x_coord, unsigned int* y_coord);
/**
 * Final steps of the scalar multiplication
 * 		(recover affine coordinates, denormalize, etc)
 * - output: output (result of the scalar multiplication)
 * - modulus_p
 * - data_a: domain parameter
 * - data_b: domain parameter
 */
void ecdh_scalar_multiplication_final(unsigned int* output, unsigned int* modulus_p, unsigned int* data_a, unsigned int* data_b);
/**
 * Turn affine into Jacobian coordinates with fixed Z
 * 		outputs in (EXPX,EXPY,EXPZ,EXPW)
 * - data_a: domain parameter
 * - x_coord: X-affine coordinate
 * - y_coord: Y-affine coordinate
 */
void affine_to_jacobian_coordinates_norand(unsigned int* data_a, unsigned int* x_coord, unsigned int* y_coord);
/**
 * Randomization of Jacobian coordinates
 * - x_coord: X-coordinate
 * - y_coord: Y-coordinate
 * - z_coord: random value (as Z-coordinate)
 */
void randomize_jacobian_coordinates(unsigned int* x_coord, unsigned int* y_coord, unsigned int* z_coord);
/**
 * Randomization of Jacobian coordinates
 * - x_coord1: X-coordinate (x_coord1 * random_s^2)
 * - x_coord2: X-coordinate (x_coord2 * random_s^2)
 * - y_coord1: Y-coordinate (y_coord1 * random_s^3)
 * - y_coord2: Y-coordinate (y_coord2 * random_s^3)
 * - ym_coord: Y-coordinate (ym_coord * random_s^3)
 * - random_s: random value (as Z-coordinate)
 */
void randomize_jacobian_segments_5(int x_coord1, int y_coord1, int x_coord2, int y_coord2, int ym_coord, unsigned int* random_s);
/**
 * Point tripling from affine in standard representation to co-Z coordinates
 *    ouputs co-Z coordinates: 
 * 		- [3]P: (1,2)
 * 		- P: (3,4)		 
 * - data_a: domain parameter
 * - x_coord: X-affine coordinate
 * - y_coord: Y-affine coordinate
 */
void coz_initial_triple(unsigned int* data_a, unsigned int* x_coord, unsigned int* y_coord);
/**
 * Conditional copy
 * 		if select_bit: dst=src1
 * 		otherwise: dst=src2
 * - dst: destination array
 * - src1: first source array
 * - src2: second source array
 */
void conditional_copy(unsigned int* dst, unsigned int *src1, unsigned int *src2, unsigned int select_bit);
/**
 * Affine from co-Z coordinates
 * 		other inputs: 
 * 			- co-Z coordinates of P in (TMP3, TMP4)
 * 			- co-Z coordinates of Q in (TMP1, TMP2)
 * 		outputs (affine coordinates of Q) in (ACCX, ACCY)
 * - data_a: domain parameter
 * - x_coord: X-affine coordinate of P
 * - y_coord: Y-affine coordinate of P
 */
void coz_coord_recovery(unsigned int* data_a, unsigned int* x_coord, unsigned int* y_coord);
/**
 * Doubling and addition 2P+Q
 * 		inputs:
 * 			- co-Z coordinates of P in (1,2)
 * 			- co-Z coordinates of Q in (3,4)
 * 		outputs:
 * 			- co-Z coordinates of 2P+Q in (1,2)
 * 			- co-Z coordinates of Q in (3,4)
 */
void coz_double_add();
/**
 * Select i-th bit of k
 * - k: array of integers
 * - i: index of the bit to select
 */
unsigned int getbit_gen(const unsigned int* k, const int i);
/**
 * Point addition
 * 		inputs:
 * 			- co-Z coordinates of P in (TMP3,TMP4)
 * 			- co-Z coordinates of Q in (TMP1,TMP2)
 * 		outputs:
 * 			- co-Z coordinates of P+Q in (TMP5,TMP6)
 * 			- co-Z coordinates of P in (ACCX,ACCY)
 */
void coz_add(void);
/**
 *  Jacobian addition
 *  -> If aux=0, then perform a Jacobian addition of the exponentiator to the i-th accumulator
 *  -> If aux=1, then perform a Jacobian addition of the j-th accumulator to the i-th accumulator
 *  -> If aux=2, then perform a Jacobian addition of the point at address in2 to the i-th accumulator
 * inputs in: in1 - index of the accumulator (i)
 *        in2 - index of the accumulator (j) / address of the point to add
 *        aux - mode in {0,1,2}
 */
void jacobian_addition(void);
/**
 *  Copy of integer arrays
 *  - in1: register to copy
 *  - out: destination register
 *  - aux: length of in1
 */
void reg_copy(void);
/**
 *  Exclusive-or of integer arrays
 *  - in1: input register
 *  - in2: input register
 *  - out: destination register
 *  - aux: length of input arrays
 */
void reg_xor(void);
/**
 *  Copy ACCX, ACCY and ACCZ in out
 *  - ACCX: input register
 *  - ACCY: input register
 *  - ACCZ: destination register
 *  - out: output register
 */
void copy_acc_reg_3(void);
/**
 * Compute the coordinate inverse
 * 		inputs in:
 * 			- in1: input to inverse
 * 			- in2: length
 * 			- aux: 
 * 		output in in1
 */
void coord_invert(void);
/**
 * Compute the Montgomery multiplication of 
 * in1 and in2 of length LEN in out
 * - in1: first source register
 * - in2: second source register
 * - out: destination register
 */
void montgomery_mult(void);
/**
 * Compute the Montgomery multiplication of 
 * in1 and in2 of length LEN in out
 * without setting the base register
 * - in1: first source register
 * - in2: second source register
 * - out: destination register
 */
void montgomery_mult_after(void);
/**
 * Compute the squaring (Montgomery multiplication)
 * of in1 of length LEN in out
 * - in1: source register
 * - out: destination register
 */
void montgomery_square(void);
/**
 * Compute the squaring (Montgomery multiplication)
 * of in1 of length LEN in out
 * without setting the base register
 * - in1: source register
 * - out: destination register
 */
void montgomery_square_after(void);
/**
 * Compute the arithmetic subtraction modulo MOD (in1-in2)
 * - in1: first source register
 * - in2: second source register
 * - out: destination register
 */
void mod_sub(void);
/**
 * Compute the modular reduction modulo MOD
 * - in1: source register
 * - out: destination register
 */
void mod_red(void);
/**
 * Recover the standard representation
 * - in1: source register
 * - out: destination register
 */
void standard_repr(void);
/**
 * Compute the Montgomery representation
 * - 3: index of a temporary register
 * - in1: source register
 * - out: destination register
 */
void montgomery_repr(void);
/**
 * Compute the arithmetic operation in2+(aux*in1)
 * - in2: first source register
 * - in1: third source register
 * - aux: index of the third source register
 * - 1: length of aux
 * - out: index of the destination register
 */
void int_macc_32(void);
/**
 * Compute the arithmetic operation in2+(aux*in1)
 * - in2: first source register
 * - in1: third source register
 * - aux: index of the third source register
 * - 1: length of aux and in2
 * - out: index of the destination register
 */
void int_macc(void);
/**
 * Compute the multiplication between in1 and in2
 * - in1: first source register
 * - in2: second source register
 * - 1: length of the inputs
 * - out: destination register
 */
void int_mult (void);
/**
 *  Compute the opposite or the identity: 
 * 			dst <- (-1)^sign * src
 *  - dst: output register
 *  - sign: sign to determine the opposite or the identity
 *  - src: source register
 */
void neg_mod_p(unsigned int* dst, unsigned char sign, unsigned int* src);
/**
 *  Copy in1 on out+2 (with out[0]=out[1]=0)
 */
void mod_import(void);
/**
 *  Copy in1+1 on out on size min(in1[0],LEN + 1)
 */
unsigned int op_import(void);
/**
 *  Copy in1 on out on size LEN
 */
void op_export(void);
/**
 *  Get len random unsigned integers in out
 * 	- len: number of unsigned integers to copy in out
 */
void random_op(unsigned int len);
/**
 *  Get a specific bit in in1
 * - in1: input register
 * - in2: index for bit selection
 */
unsigned char get_bit(void);
/**
 *  Get a specific digit from in1
 * - in1: input register
 * - in2: index for digit selection
 */
unsigned char get_digit_l2r(void);
/**
 * Load the modulus in the coprocessor
 * - in1+2: modulus
 * - aux: length of the modulus
 */
void set_mod(void);
/**
 * Return a random unsigned integer from the TRNG 
 * or from the external source, according to the
 * value of ECDH_RANDOM_SOURCE
 */
unsigned int ecdh_get_random(void);
/**
 * Euclidean division
 * - scalar: initial scalar
 * - div_32: 64-bit random
 * - quo: division quotient
 * - rem: division reminder
 * - div: updated random
 */
void eucl_div_64(unsigned int* scalar, unsigned int* div_32, unsigned int* quo, unsigned int* rem, unsigned int* div);
/**
 * Euclidean division scalar/div
 * - quo: quotient
 * - op: scalar
 * - div: to divide with
 */
unsigned int eucl_div_32(unsigned int* quo, unsigned int* op, unsigned int div);

/******************************************************************************/
/***                             CONSTANTS                                  ***/
/******************************************************************************/

/** mask for conditional copies */
unsigned int mask[2] = {0x66666666, 0x99999999};

/******************************************************************************/
/***                             RANDOMNESS                                 ***/
/******************************************************************************/

unsigned int ecdh_get_random()
{		
	unsigned int rand;
	switch (ECDH_RANDOM_SOURCE)
	{
		case ECDH_INTERNAL_RANDOM_SOURCE:
			rand = random_word() & 0xFFFFFFFF;
			break;
			
		case ECDH_EXTERNAL_RANDOM_SOURCE:
			rand = ECDH_RANDOM_BUFFER[delta] & 0xFFFFFFFF;
			delta++;
			break;
	}
	return rand;
}

/******************************************************************************/
/***                           ECDH primitive                               ***/
/******************************************************************************/

void ecdh(int level, unsigned int* output, unsigned int* modulus, unsigned int* data_a, unsigned int* data_b, unsigned int* data_r, unsigned int* cm_param, unsigned int* scalar, unsigned int* x_coord, unsigned int* y_coord, int len)
{
	unsigned char param;

	//--- Switch on CCP ---//
	copro_on();

	//--- Switch on TRNG ---//
	if (ECDH_RANDOM_SOURCE == ECDH_INTERNAL_RANDOM_SOURCE)
	{
		RNG_on();
	} 
	else
	{
		delta =0;
	}

	nb_reg = (unsigned char)(*cm_param);	
	param = nb_reg;
	if (nb_reg < 2)
		nb_reg =2;

	//--- Scalar multiplication ---//

	ecdh_scalar_multiplication_init(modulus,len);
	ecdh_scalar_multiplication_core(level, param, data_a, data_r, scalar, x_coord, y_coord);
	ecdh_scalar_multiplication_final(output, modulus, data_a, data_b);

	//--- Switch off TRNG ---//
	if (ECDH_RANDOM_SOURCE == ECDH_INTERNAL_RANDOM_SOURCE)
	{
		RNG_off();
	}

	// -- Switch off CCP --//
	copro_off();
}

/******************************************************************************/
/***                           ECDH main functions                          ***/
/******************************************************************************/

void ecdh_scalar_multiplication_init(unsigned int* modulus, int len)
{
	//--- Load modulus ---//
	in1 = modulus;
	out = MOD;
	mod_import ();

	//--- Montgomery precomputation (computation of 1/K^2 mod MOD, K = 2^{-j} mod m, with j = 32*6) ---//
	
	// Compute -p^(-1) mod 2^32
	in1 = MOD;
	aux = len;
	set_mod ();

	// Compute 1/K^2 mod MOD
	low_montgomery_inv(1,2);
}

void ecdh_scalar_multiplication_core(int level, unsigned char param, unsigned int* data_a, unsigned int* data_r, unsigned int* scalar, unsigned int* x_coord, unsigned int* y_coord)
{
	int b,i;

	if (level==1)
	{
		// CM order mult
		SCAL1[0] = ecdh_get_random();
		b = scalar[0] & 1; // b <- k mod 2

		SCAL[LEN]=0;
		in1 = data_r+1;
		in2 = scalar;
		aux = SCAL1[0];
		out = SCAL; 
		int_macc_32();

		ecdh_scalar_point_mult_binary_coz(level,param,data_a,data_r,SCAL,LEN+1,x_coord,y_coord);

		affine_to_jacobian_coordinates_norand(data_a,ACCX,ACCY);
   
	    aux = LEN;
	    in1 = EXPX;
	    out = ACCX;
	    reg_copy();	

	    in1 = EXPY;
	    out = ACCY;
	    reg_copy();

		in1 = EXPZ;
	    out = ACCZ;
	    reg_copy(); 
	} 
	else
	{
		// case level=2
		// Generate a 64-bit random
		SCAL[0]=0x80000000 | ecdh_get_random();
		SCAL[1]=0x80000000 | ecdh_get_random();
		

		//  -- Euclidean Division
		eucl_div_64(scalar, SCAL, SCAL3, SCAL1, SCAL2);
		// SCAL3 = floor(k/random), SCAL1 = k mod random, SCAL2 = random

		//	-- Compute [random].P
		ecdh_scalar_point_mult_binary_coz(level,param,data_a,data_r,SCAL2,2,x_coord,y_coord);

		//	-- Compute [k/random].([random].P)

		// affine coordinates in ACCX,ACCY
		aux=LEN;
		in1 = ACCX;
		out = AFF_X;
		reg_copy();
		in1 = ACCY;
		out = AFF_Y;
		reg_copy();

		ecdh_scalar_point_mult_binary_coz(level,param,data_a,data_r,SCAL3,LEN,AFF_X,AFF_Y);

		affine_to_jacobian_coordinates_norand(data_a,ACCX,ACCY);

		aux = LEN;
	    in1 = EXPX;
	    out = JAC_X;
	    reg_copy();	

	    in1 = EXPY;
	    out = JAC_Y;
	    reg_copy();

		in1 = EXPZ;
	    out = JAC_Z;
	    reg_copy(); 

		//	-- Compute [k mod random].P

		ecdh_scalar_point_mult_binary_coz(level,param,data_a,data_r,SCAL1,2,x_coord,y_coord);

		affine_to_jacobian_coordinates_norand(data_a,ACCX,ACCY);

	    aux = LEN;
	    in1 = EXPX;
	    out = ACCX;
	    reg_copy();	

	    in1 = EXPY;
	    out = ACCY;
	    reg_copy();

		in1 = EXPZ;
	    out = ACCZ;
	    reg_copy(); 

		//	-- Add [k/random].([random].P) and [k mod random].P
		aux = 2;
		in1 = 0;
		in2 = JAC_X;
		jacobian_addition();
	}
}

void ecdh_scalar_multiplication_final(unsigned int* output, unsigned int* modulus_p, unsigned int* data_a, unsigned int* data_b)
{
	unsigned char i, j;
	
	// Copy result at begining of workzone after the modulus
	out = TMPX;
	COPY_ACC_LEN_3 ();
	
	//--- Invert Z-coord ---//
	in1 = TMPZ;
	in2 = ((unsigned int*)LEN);
	aux = 3 * NBACC + 7;
	INVERT (); 

	//--- Recover affine coordinates ---//

	// Compute Z^(-2) and store it in TMP2
	
	in1 = TMPZ;
	out = TMP2;
	
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(4,2,2);
	low_get_seg(4,out,LEN);

	// Compute Z^(-3) and store it in TMP1
	
	in1 = TMPZ;
	in2 = TMP2;
	out = TMP1;
	
	montgomery_mult ();

	in1 = TMPY;
	in2 = TMP1;
	out = RESY;
	
	montgomery_mult ();
	
	in1 = TMPX;
	in2 = TMP2;
	out = RESX;
	
	montgomery_mult ();

	//--- Compute (y^2 - x^3 - a.x) and load b for FA check ---//

	// reload modulus for FA checks 	

	in1 = modulus_p;
	out = MOD;
	mod_import ();
	
	in1 = MOD;
	aux = LEN;
	set_mod ();

	// TMP1 = y^2/K
	
	in1 = RESY;
	out = TMP1;
	
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(4,2,2);
	low_get_seg(4,out,LEN);

	// TMP2 = x^2/K
	
	in1 = RESX;
	out = TMP2;
	
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(4,2,2);
	low_get_seg(4,out,LEN);
	
	// TMP3 = x^3/K
	
	in1 = RESX;
	in2 = TMP2;
	out = TMP3;
	
	montgomery_mult ();
	
	// TMP2 = (y^2 - x^3)/K
	
	in1 = TMP1;
	in2 = TMP3;
	out = TMP2;
	
	mod_sub (); 	

	// Denormalize: TMP1 = y^2 - x^3
	
	in1 = TMP2;
	out = TMP1;
	
	standard_repr (); 
	
	// TMP3 = a.x
	in1 = data_a;
	out = TMP2;
	
	op_import ();
	
	in1 = TMP2;
	in2 = RESX;
	out = TMP3;

	montgomery_mult (); 

	// TMP2 = y^2 - x^3 - a.x
	
	in1 = TMP1;
	in2 = TMP3;
	out = TMP2;
	
	mod_sub ();


	// complete reduction
	
	in1 = TMP2;
	out = TMP1;
	mod_red ();

	// TMP1 = b
	
	in1 = data_b;
	out = TMP2;
	
	op_import ();
	
	//--- Check, denormalize and output point ---//

	// 1st check : if (y^2 - x^3 - a.x != b) erase x- and y-coordinates

	for (i = 0; i < OPSIZE(LEN); i++) 
	{
		if (((TMP1[i]^TMP2[i])&0xFFFFFFFF) != 0)
		{
			// erase x-coordinate
			for (j = 0; j < OPSIZE(LEN); j++) *(RESX + j) = 0;
			
			// erase y-coordinate
			for (j = 0; j < OPSIZE(LEN); j++) *(RESY + j) = 0;
		}
	}
		
	// denormalize and output point
	// X-coord
	// denormalize

	in1 = RESX;
	out = TMP3;
	
	standard_repr ();

	// output
	
	in1 = TMP3;
	out = RESX;
	mod_red ();

	in1 = RESX;
	out = output+LEN;
	
	op_export ();
	
	// Y-coord
	// denormalize
	
	in1 = RESY;
	out = TMP4;
	
	standard_repr ();
	
	// output
	
	in1 = TMP4;
	out = RESY;
	mod_red ();

	in1 = RESY;
	out = output;
	
	op_export ();
	
	// 2nd check : if (y^2 - x^3 - a.x != b) erase output

	for (i = 0; i < OPSIZE(LEN); i++) 
	{
		if (((TMP1[i]^TMP2[i])&0xFFFFFFFF) != 0)
		{
			// erase x-coordinate
			for (j = 0; j < LEN; j++) output[j] = 0;
			
			// erase y-coordinate
			for (j = 0; j < LEN; j++) output[j+LEN] = 0;
		}
	}
	
}

/******************************************************************************/
/***                          Scalar Multiplications                        ***/
/******************************************************************************/

/***                          Binary COZ 						***/

// swap 4 and 7 segments if b=1 using 5
void swap (unsigned int b)
{
	low_copy_seg(5,7);
	low_copy_seg(4+3*(b&1),4);
	low_copy_seg(7-3*(b&1),5);
}

void ecdh_scalar_point_mult_binary_coz(unsigned char param, unsigned int* data_a, unsigned int* data_r, unsigned int *scalar, unsigned int scalar_len, unsigned int* x_coord, unsigned int* y_coord)
{    
    int i,j,n;
    unsigned int b;
    n = scalar_len*32;
    
    // set n to msb index    
    while (getbit_gen(scalar, n-1) != 1) { n--; }

    in1= scalar;
    aux = scalar_len;
    out = SCAL;
    reg_copy();

#ifdef ECC_CM
    for (j=1; j<param+1; j++)
	{		
		// countermeasure refresh & swap
		out = S_MASK+(j-1)*(LEN+1);
    	random_op(LEN);

    	in1= SCAL;
	    in2= S_MASK+(j-1)*(LEN+1);
	    aux = scalar_len;
	    out = SCAL;
	    reg_xor();
	}
#endif
  
    coz_initial_triple(data_a,x_coord,y_coord);

    low_get_seg(1,TMP1,LEN);	// R1x
    low_get_seg(2,TMP2,LEN);	// R1y
    low_get_seg(3,TMP3,LEN);	// R0x
    low_get_seg(4,TMP4,LEN);	// R0y

    out = TMP7;
    random_op(LEN);
    randomize_jacobian_coordinates(TMP1,TMP2,TMP7);
    randomize_jacobian_coordinates(TMP3,TMP4,TMP7);

    low_set_seg(1,TMP1,LEN);
    low_set_seg(2,TMP2,LEN);
    low_set_seg(3,TMP3,LEN);
    low_set_seg(4,TMP4,LEN);

    // main loop
    for (i=n-2; i>=1; i--)
    {
   		low_mod_sub_mod(7,4);	

#ifdef ECC_CM				
		out = TMP7;
		random_op(LEN);
		randomize_jacobian_segments_5(1, 2, 3, 4, 7, TMP7);			
#endif	
			
   		b= getbit_gen(SCAL, i) ^ getbit_gen(SCAL, i+1);
   		swap(b);

#ifdef ECC_CM		
		for (j=1; j<param+1; j++)
		{
			// refresh
				
			out = TMP7;
			random_op(LEN);
			randomize_jacobian_segments_5(1, 2, 3, 4, 7, TMP7);

			b=getbit_gen(S_MASK+(j-1)*(LEN+1), i) ^ getbit_gen(S_MASK+(j-1)*(LEN+1), i+1);
			swap(b);
		}
#endif

      // co-Z double and add
      coz_double_add();	

    }
    low_get_seg(1,TMP1,LEN);
	low_get_seg(2,TMP2,LEN);
	low_get_seg(3,TMP3,LEN);
	low_get_seg(4,TMP4,LEN);

    // final addition

    b = getbit_gen(SCAL, 1);    
#ifdef ECC_CM
	for (j=1; j<param+1; j++)
	{
		b ^= getbit_gen(S_MASK+(j-1)*(LEN+1), 1);
	}
#endif
		
    neg_mod_p(TMP4, b, TMP4); 	
    coz_add();
  
    // conditional copy
    b = getbit_gen(SCAL, 0);
#ifdef ECC_CM
    for (j=1; j<param+1; j++)
		{
			b ^= getbit_gen(S_MASK+(j-1)*(LEN+1), 0);
		}
#endif

    // result from coz_add or not depending on b
    conditional_copy(TMP1, TMP5, TMP1, b);
    conditional_copy(TMP2, TMP6, TMP2, b);
    conditional_copy(TMP3, ACCX, TMP3, b);
    conditional_copy(TMP4, ACCY, TMP4, b);
            
    // correct sign of Q
    neg_mod_p(TMP4, 1, TMP4);		
  
    coz_coord_recovery(data_a,x_coord,y_coord);
}

void coz_initial_triple(unsigned int* data_a, unsigned int* x_coord, unsigned int* y_coord)
{
	
	in1 = x_coord;
	low_set_seg_base(2,in1,LEN);
	low_montgomery_repr (1,2,3);

	in1 = y_coord;
	low_set_seg_base(4,in1,LEN);
	low_montgomery_repr (2,4,3);

	//--- Initial co-Z doubling ---//

	// X1^2
	low_montgomery_mult(3,1,1);

	// 2 X1^2
	low_mod_doubl(4,3);

	// 3 X1^2
	low_mod_add(7,3,4);
	low_copy_seg(3,7);

	// 3 X1^2 + a = B
	in1 = data_a+1;
	low_set_seg_base(5,in1,LEN);
	low_montgomery_repr (4,5,7);

	low_mod_add(7,3,4);
	low_copy_seg(3,7);

	// Y1^2
	low_montgomery_mult(4,2,2);

	// 2 Y1^2
	low_mod_doubl(7,4);
	low_copy_seg(4,7);

	// 4 Y1^2
	low_mod_doubl(5,4);

	// 4 X1 Y1^2 = A = X1'		
	low_montgomery_mult(7,5,1);	
	low_copy_seg(5,7);	

	// B^2
	low_montgomery_mult(6,3,3);	

	// B^2-A
	low_mod_sub(7,6,5);
	low_copy_seg(6,7);	

	// X2
	low_mod_sub(7,6,5);
	low_copy_seg(6,7);

	// A-X2
	low_mod_sub(1,5,6);

	// B(A-X2)
	low_montgomery_mult(7,1,3);
	low_copy_seg(1,7);

	// 4 Y1^4	
	low_montgomery_mult(7,4,4);

	// 8 y1^4 = Y1'
	low_mod_doubl(3,7);

	// Y2
	low_mod_sub(7,1,3);
	low_copy_seg(1,7);

	// [2]P = (T6,T1)
	// P = (T5,T3)

	low_copy_seg(4,3);
	low_copy_seg(3,5);
	low_copy_seg(2,1);
	low_copy_seg(5,6);

	// [2]P = (T5,T2)
	// P = (T3,T4)

    //--- Co-Z addition ---//
    
    // L1. X2 - X1 
    low_mod_sub(1,5,3);
    
    // L1. (X2 - X1)^2 
    low_montgomery_mult(7,1,1);
    low_copy_seg(1,7);
    
    // L2. X1' = X1 (X2 - X1)^2 
    low_montgomery_mult(7,1,3);
    low_copy_seg(3,7);
    
    // L3. X2 (X2 - X1)^2
    low_montgomery_mult(7,5,1);
    low_copy_seg(5,7);
    
    // L4. Y2 - Y1
    low_mod_sub(7,2,4);
    low_copy_seg(2,7);
    
    // L4. (Y2 - Y1)^2
    low_montgomery_mult(1,2,2);
    
    // L6. (Y2 - Y1)^2 - X1 (X2 - X1)^2
    low_mod_sub(7,1,3);
    low_copy_seg(1,7);

    // L6. X3 = (Y2 - Y1)^2 - (X1 + X2) (X2 - X1)^2
    low_mod_sub(7,1,5);
    low_copy_seg(1,7);
    
    // L5. C-B= (X2 - X1)^3
    low_mod_sub(7,5,3);
    low_copy_seg(5,7);
     
    // L5. Y1' = Y1 (X2 - X1)^3
    low_montgomery_mult(7,4,5);
    low_copy_seg(4,7);
      
    // X1 (X2 - X1)^2 - X3
    low_mod_sub(5,3,1);
    
    // (Y2 - Y1) (X1 (X2 - X1)^2 - X3)
    low_montgomery_mult(7,2,5);
    
    // Y3 = (Y2 - Y1) (X1 (X2 - X1)^2 - X3) - Y1'
    low_mod_sub(2,7,4);

    // (X3,Y3) = (TMP1,TMP2)
    // (X1',Y1') = (TMP3,TMP4)
}

void coz_double_add()
{
 
    //--- first addition ---//
    
    // X2 - X1
	low_mod_sub(5,3,1);	

    // (X2 - X1)^2
	low_montgomery_mult_after(7,5,5);
	low_copy_seg(5,7);
    
    // X1' = X1 (X2 - X1)^2
	low_montgomery_mult_after(7,1,5);
	low_copy_seg(1,7);
    
    // X2 (X2 - X1)^2
	low_montgomery_mult_after(7,3,5);
	low_copy_seg(3,7);

    // Y2 - Y1
	low_mod_sub(7,4,2);	
	low_copy_seg(4,7);
    
    // (Y2 - Y1)^2
	low_montgomery_mult_after(5,4,4);
    
    // (Y2 - Y1)^2 - X1 (X2 - X1)^2
	low_mod_sub(6,5,1);	
    
    // X3 = (Y2 - Y1)^2 - (X1 + X2) (X2 - X1)^2
	low_mod_sub(7,6,3);	
	low_copy_seg(6,7);
    
    // Y2 - Y1 + X1'
	low_mod_add(7,4,1);	
	low_copy_seg(4,7);		
    
    // Y2 - Y1 + X1' - X3
	low_mod_sub(7,4,6);	
    
    // (Y2 - Y1 + X1' - X3)^2
	low_montgomery_mult_after(4,7,7);
    
    // (Y2 - Y1 + X1' - X3)^2 - (Y2 - Y1)^2
	low_mod_sub(7,4,5);	
	low_copy_seg(4,7);
    
    // X1' - X3
	low_mod_sub(5,1,6);	
    
    // (X1' - X3)^2
	low_montgomery_mult_after(7,5,5);
	low_copy_seg(5,7);
    
    // 2 (Y2 - Y1) (X1' - X3)^2
	low_mod_sub(7,4,5);	
	low_copy_seg(4,7);
    
    // (X2 - X1)^3
	low_mod_sub(7,3,1);	
	low_copy_seg(3,7);
    
    // Y1 (X2 - X1)^3
	low_montgomery_mult_after(7,2,3);
    
    // 2 Y1' = 2 Y1 (X2 - X1)^3
	low_mod_add(2,7,7);			
    
    // 2 Y3 = (Y2 - Y1) (X1' - X3)^2 - 2 Y1 (X2 - X1)^3
	low_mod_sub(7,4,2);	
	low_copy_seg(4,7);
    
    //--- conjugate addition ---//
    
    // X3 (X1' - X3)^2
	low_montgomery_mult_after(7,5,6);
	low_copy_seg(6,7);
    
    // 4 X3 (X1' - X3)^2
	low_mod_doubl(7,6);
	low_mod_doubl(6,7);
    
    // X1' (X1' - X3)^2
	low_montgomery_mult_after(7,1,5);	
	low_copy_seg(1,7);
    
    // 4 X1' (X1' - X3)^2
	low_mod_doubl(7,1);
	low_mod_doubl(5,7);
    
    // 2 (Y1' - Y3)
	low_mod_sub(1,2,4);	
    
    // 4 (Y1' - Y3)^2
	low_montgomery_mult_after(7,1,1);
    
    // 4 (Y1' - Y3)^2 - 4 X3 (X1' - X3)^2
	low_mod_sub(1,7,6);	
    
    // 4 X4 = 4 (Y1' - Y3)^2 - 4 (X1' + X3) (X1' - X3)^2
	low_mod_sub(7,1,5);	
	low_copy_seg(1,7);
    
    // 2 (Y1' + Y3)
	low_mod_add(3,2,4);			
    
    // 4 (Y3 + Y1')^2
	low_montgomery_mult_after(7,3,3);
    
    // 4 (Y3 + Y1')^2 - 4 X3 (X1' - X3)^2
	low_mod_sub(3,7,6);	
    
    // 4 X2' = 4 (Y3 + Y1')^2 - 4 (X1' + X3) (X1' - X3)^2
	low_mod_sub(7,3,5);	
	low_copy_seg(3,7);
    
    // 4 (X1' - X3)^3
	low_mod_sub(7,5,6);	
    
    // 8 Y3 (X1' - X3)^3
	low_montgomery_mult_after(5,4,7);

    // 2 (Y1' - Y3)
	low_mod_sub(7,2,4);	
	low_copy_seg(2,7);
    
    // 4 Y3
	low_mod_doubl(7,4);
    
    // 2 (Y1' + Y3)
	low_mod_add(4,7,2);			
    
    // 4 X2' - 4 X3 (X1' - X3)^2
	low_mod_sub(7,3,6);	
	low_copy_seg(3,7);
    
    // 8 (Y1' + Y3) [X2' - X3 (X1' - X3)^2)
	low_montgomery_mult_after(7,3,4);
    
    // 8 Y2' = 8 (Y1' + Y3) [X2' - X3 (X1' - X3)^2) - 8 Y3 (X1' - X3)^3
	low_mod_sub(4,7,5);	
    
    // 4 X2' 
	low_mod_add(7,3,6);		
	low_copy_seg(3,7);	
    
    // 4 [X4 - X3 (X1' - X3)^2) 
	low_mod_sub(7,6,1);	
	low_copy_seg(6,7);	
    
    // 8 (Y1' - Y3) [X4 - X3 (X1' - X3)^2) 
	low_montgomery_mult_after(7,2,6);
    
    // 8 Y4 = 8 (Y1' - Y3) [X4 - X3 (X1' - X3)^2) - 8 Y3 (X1' - X3)^3
	low_mod_sub(2,7,5);	
}

unsigned int getbit_gen(const unsigned int* k, const int i)
{
    unsigned int W,b;
    
    unsigned int qi = i/32;
    unsigned int ri = i%32;
    
    W = k[qi];
    b = (W >> ri) & 1;
    
    return b;
}

void conditional_copy(unsigned int* dst, unsigned int *src1, unsigned int *src2, unsigned int select_bit)
{
    int i;
    unsigned int tmp;
    
    for(i=0; i<LEN; i++)
    {
        tmp = (src1[i] & 0x99999999) ^ (src2[i] & 0x66666666);
        tmp = tmp ^ (src1[i] & mask[select_bit]);
        tmp = tmp ^ (src2[i] & mask[select_bit]);
        
        dst[i] = tmp;
    }
}

void neg_mod_p(unsigned int* dst, unsigned char sign, unsigned int* src)
{
	// dst <- (-1)^sign * src
	in1 = src;
	out = TMP7;
	low_set_seg_base(2,in1,LEN);
	low_mod_sub_mod(4,2);	
	low_get_seg(4,out,LEN);
	if ((sign&1)==0)
		in1 = src;
	else
		in1 = TMP7;
	out = dst;
	aux = LEN;
	reg_copy();
	
}

void coz_add()
{
	in1 = TMP1;
    in2 = TMP3;
    out = TMP5;
    mod_sub();
    in1 = TMP5;
    montgomery_square();
    in1 = TMP3;
    in2 = TMP5;
    out = ACCX;
    montgomery_mult();
    in1 = TMP1;
    out = ACCZ;
    montgomery_mult();
    in1 = TMP2;
    in2 = TMP4;
    out = TMP6;
    mod_sub();
    in1 = TMP6;
    out = TMP5;
    montgomery_square_after();
    in1 = TMP5;
    in2 = ACCX;
    out = TMP5;
    mod_sub();
    in2 = ACCZ;
    mod_sub();
    in1 = ACCZ;
    in2 = ACCX;
    out = ACCZ;
    mod_sub();
    in1 = TMP4;
    in2 = ACCZ;
    out = ACCY;
    montgomery_mult_after();
    in1 = ACCX;
    in2 = TMP5;
    out = ACCZ;
    mod_sub();
    in1 = TMP6;
    in2 = ACCZ;
    out = TMP6;
    montgomery_mult_after();
    in2 = ACCY;
    mod_sub();
}

void coz_coord_recovery(unsigned int* data_a, unsigned int* x_coord, unsigned int* y_coord)
{
	in1 = x_coord;
	out = TMP5;
	montgomery_repr();

	in1 = y_coord;
	out = SCAL;
	montgomery_repr();	

	// x0 * Y0
	in1 = TMP5;
	in2 = TMP4;
	out = TMP5;
	montgomery_mult();

	// (x0 * Y0)^-1

	in1 = TMP5;
	in2 = ((unsigned int*)LEN);
	aux = 3 * NBACC + 7;
	INVERT (); 

	in1 = TMP5;
	in2 = TMP3;
	out = TMP5;
	montgomery_mult();
	in1 = TMP5;
	in2 = SCAL;
	out = TMP5;
	montgomery_mult();
	in1 = TMP5;
	out = TMP6;
	montgomery_square();
	in1 = TMP6;
	in2 = TMP5;
	out = TMP5;
	montgomery_mult();
	in1 = TMP1;
	in2 = TMP6;
	out = ACCX;
	montgomery_mult();
	in1 = TMP2;
	in2 = TMP5;
	out = ACCY;
	montgomery_mult();

	in1 = ACCX;
	out = ACCX;
	standard_repr();

	in1 = ACCY;
	out = ACCY;
	standard_repr();
}

void randomize_jacobian_coordinates(unsigned int* x_coord, unsigned int* y_coord, unsigned int* random_s)
{
	
	low_set_seg_base(5,random_s,LEN);
	
	// S^2
	low_montgomery_mult(6,5,5);

	// X = X*S^2
	low_set_seg_base(7,x_coord,LEN);
	low_montgomery_mult(5,6,7);
	low_get_seg(5,x_coord,LEN);

	// S^3
	low_set_seg_base(5,random_s,LEN);
	low_montgomery_mult(7,5,6);

	// Y <- Y*S^3
	low_set_seg_base(6,y_coord,LEN);
	low_montgomery_mult(5,7,6);
	low_get_seg(5,y_coord,LEN);
}

void randomize_jacobian_segments_5(int x_coord1, int y_coord1, int x_coord2, int y_coord2, int ym_coord, unsigned int* random_s)
{
	
	low_set_seg_base(5,random_s,LEN);
	
	// S^2
	low_montgomery_mult(6,5,5);

	// X1 = X1*S^2
	low_montgomery_mult(5,6,x_coord1);
	low_copy_seg(x_coord1,5);

	// X2 = X2*S^2
	low_montgomery_mult(5,6,x_coord2);
	low_copy_seg(x_coord2,5);

	// S^3
	low_set_seg_base(5,random_s,LEN);
	low_montgomery_mult(6,5,6);

	// Y1 <- Y1*S^3
	low_montgomery_mult(5,6,y_coord1);
	low_copy_seg(y_coord1,5);

	// Y2 <- Y2*S^3
	low_montgomery_mult(5,6,y_coord2);
	low_copy_seg(y_coord2,5);

	// -Y <- -Y*S^3
	low_montgomery_mult(5,6,ym_coord);
	low_copy_seg(ym_coord,5);
}

void affine_to_jacobian_coordinates_norand(unsigned int* data_a, unsigned int* x_coord, unsigned int* y_coord)
{
	int i;

	// Define Z-coordinate (Z/K)
	
	out = EXPZ;
	EXPZ[0]= 1;
	for (i=1; i<LEN; i++)
	{
		EXPZ[i]=0;
	}
	low_set_seg_base(4,EXPZ,LEN);
	low_montgomery_inv(1,3);
	low_montgomery_mult(2,4,1);		// 2 = Z/K
	low_get_seg(2,EXPZ,LEN);

	// Compute Z^2/K^2 

	low_montgomery_mult(4,2,1);		// Z/K^2
	
	// Z^2/K^2
	
	low_montgomery_mult(5,2,4);		// Z^2/K^2
	
	// Compute Z^3/K^2 (stored in TMP3)
	
	low_montgomery_mult(4,2,5);

	// Load, randomize and normalize X-coordinate
	// Load X-coordinate
	
	// X: secret key X coordinate
	// Compute X.Z^2/K

	low_set_seg_base(6,x_coord,LEN);
	low_montgomery_mult(3,6,5);
	low_get_seg(3,EXPX,LEN);

	// Load, randomize and normalize Y-coordinate
	// Load Y-coordinate
	
	// Y: secret key Y coordinate (data_Gy)
	// Compute Y.Z^3/K

	low_set_seg_base(6,y_coord,LEN);
	low_montgomery_mult(5,6,4);
	low_get_seg(5,EXPY,LEN);

	// Compute normalized W-coordinate (= a.Z^4/K)
	// Load a
	// Compute a.Z^3/K

	low_set_seg_base(6,data_a+1,LEN);
	low_montgomery_mult(3,6,4);
	
	// Compute a.Z^4/K
	low_montgomery_mult(4,3,2);
	low_get_seg(4,EXPW,LEN);

	// End of Jacobian coordinates computation
	// Here (EXPX,EXPY,EXPZ,EXPW) contains Jacobian coordinates (x_coord,Y_coord)
}

/******************************************************************************/
/***                          Sub Fonctions                                 ***/
/******************************************************************************/

unsigned int eucl_div_32(unsigned int* quo, unsigned int* op, unsigned int div)
{
    int i;
    
    unsigned int q;
    u64 tmp1, tmp2;
    unsigned int tmpRemH;
    
    q = op[LEN-1]/div;
    tmpRemH = op[LEN-1] - q*div;
    quo[LEN-1] = q;
    
    for (i=LEN-2; i>=0; i--)
    {
        tmp1 = (((u64) tmpRemH) << 32) + op[i];
        q = (unsigned int) (tmp1 / div);
        tmp2 = tmp1 - ((u64) q) * ((u64) div);
        tmpRemH = (unsigned int) tmp2;
        quo[i]=q;
    }
    
    return tmpRemH;
}

void eucl_div_64(unsigned int* scalar, unsigned int* div_32, unsigned int* quo, unsigned int* rem, unsigned int* div)
{
    unsigned int rem_32;
    int i;

    rem[0]=eucl_div_32(quo, scalar, div_32[0]);
	for (i=1; i<LEN; i++)
    {
    	rem[i]=0;
    }

	rem_32=eucl_div_32(quo, quo, div_32[1]);

	// rem = rem + div*rem_32
	in1 = rem;
	in2 = div_32;
	aux = rem_32;
	out = rem;
	int_macc();

	// updated div
	in1 = div_32;
	in2 = div_32+1;
	out = div;
	int_mult();
}

void jacobian_addition()
{
	unsigned char accInd;
	unsigned char product;
	unsigned int *X2, *Y2, *Z2;
	unsigned int *X1, *Y1, *Z1;

	accInd = ((unsigned int)in1);

	product = accInd * OPSIZE(LEN);
	X1 = ACCX + product;
	Y1 = ACCY + product;
	Z1 = ACCZ + product;

	switch (aux)
	{
		case 0:
			// mode 0
			
			X2 = EXPX;
			Y2 = EXPY;
			Z2 = EXPZ;
			break;
			
		case 1:
			// mode 1
			
			X2 = ACCX + ((unsigned int)in2) * OPSIZE(LEN); 
			Y2 = ACCY + ((unsigned int)in2) * OPSIZE(LEN);
			Z2 = ACCZ + ((unsigned int)in2) * OPSIZE(LEN);

			break;
			
		case 2:
			// mode 2
			
			X2 = in2;
			Y2 = in2 + OPSIZE(LEN);
			Z2 = in2 + 2 * OPSIZE(LEN);
			break;
	}

	// TMP4=Z1^2
	in1 = Z1;
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(4,2,2);

	// TMP1=X2*TMP4

	out = TMP1;
	in1 = X2;
	low_set_seg_base(5,in1,LEN);
	low_montgomery_mult(3,5,4);
	low_get_seg(3,out,LEN);

	// TMP5=Y2*Z1

	in1 = Y2;

	low_set_seg_base(6,in1,LEN);
	low_montgomery_mult(3,2,6);

	// TMP2=TMP5*TMP4

	out = TMP2;
	low_montgomery_mult(5,3,4);
	low_get_seg(5,out,LEN);

	// TMP4=Z2^2
	in1 = Z2;

	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(4,2,2);

	// TMP5=X1*TMP4

	out = TMP5;
	in1 = X1;

	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(5,2,4);
	low_get_seg(5,out,LEN);
	
	// X1=Y1*Z2

	out = X1;
	in1 = Y1;
	in2 = Z2;

	low_set_seg_base(2,in1,LEN);
	low_set_seg_base(3,in2,LEN);
	low_montgomery_mult(6,2,3);
	low_get_seg(6,out,LEN);

	// Y1=X1*TMP4

	out = Y1;
	low_montgomery_mult(3,6,4);
	low_get_seg(3,out,LEN);

	// X1=TMP5-TMP1

	out = X1;
	in1 = TMP5;
	in2 = TMP1;

	mod_sub ();

	in1 = X1;
	out = EXPW;

	mod_red ();

	// TMP5=Z2*Z1
	in1 = Z2;
	in2 = Z1;

	low_set_seg_base(2,in1,LEN);
	low_set_seg_base(3,in2,LEN);
	low_montgomery_mult(4,2,3);

	// Z1=X1*TMP5

	out = Z1;
	in1 = X1;
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(3,2,4);
	low_get_seg(3,out,LEN);

	// TMP5=Y1-TMP2
	in1 = Y1;
	in2 = TMP2;
	low_set_seg_base(6,in1,LEN);
	low_set_seg_base(4,in2,LEN);
	low_mod_sub(3,6,4);	

	// TMP4=X1^2
	low_montgomery_mult(4,2,2);

	// TMP3=TMP5^2
	low_montgomery_mult(6,3,3);

	// Y1=TMP1*TMP4

	out = Y1;
	in1 = TMP1;
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(5,2,4);
	low_get_seg(5,out,LEN);

	// TMP1=TMP4*X1
	in2 = X1;
	low_set_seg_base(2,in2,LEN);
	low_montgomery_mult(7,4,2);

	// TMP3=TMP3-TMP1
	low_mod_sub(4,6,7);

	// TMP4=Y1+Y1
	in1 = Y1;
	low_set_seg_base(2,in1,LEN);
	low_mod_doubl(6,2);

	// X1=TMP3-TMP4

	out = X1;
	low_mod_sub(5,4,6);
	low_get_seg(5,out,LEN);

	// TMP3=Y1-X1
	low_mod_sub(6,2,5);


	// TMP6=TMP5*TMP3
	low_montgomery_mult(5,3,6);

	// TMP4=TMP2*TMP1
	in1 = TMP2;
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(6,2,7);

	// Y1=TMP6-TMP4

	out = Y1;
	low_mod_sub(4,5,6);
	low_get_seg(4,out,LEN);

}

void montgomery_mult()
{
	low_set_seg_base(2,in1,LEN);
	low_set_seg_base(3,in2,LEN);
	low_montgomery_mult(4,2,3);
	low_get_seg(4,out,LEN);
}

void montgomery_mult_after()
{
	low_set_seg_base(2,in1,LEN);
	low_set_seg_base(3,in2,LEN);
	low_montgomery_mult_after(4,2,3);
	low_get_seg(4,out,LEN);
}

void montgomery_square()
{
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult(4,2,2);
	low_get_seg(4,out,LEN);
}

void montgomery_square_after()
{
	low_set_seg_base(2,in1,LEN);
	low_montgomery_mult_after(4,2,2);
	low_get_seg(4,out,LEN);
}

void mod_sub()
{
	low_set_seg_base(2,in1,LEN);
	low_set_seg_base(3,in2,LEN);
	low_mod_sub(4,2,3);	
	low_get_seg(4,out,LEN);
}

void int_macc_32 ()	
{
	low_set_seg_base(2, in1, LEN);
	low_set_seg_base(3, in2, LEN);
	low_set_seg(4, &aux, 1);				
	low_int_macc(5, 3, LEN, 4, 1, 2, LEN);
	low_get_seg(5, out, LEN+1);		
}

void int_macc ()	
{
	low_set_seg_base(2, in1, LEN);
	low_set_seg_base(3, in2, 1);
	low_set_seg(4, &aux, 1);				
	low_int_macc(5, 2, LEN, 4, 1, 3, 1);
	low_get_seg(5, out, LEN);		
}

void int_mult ()
{
	low_set_seg_base(2,in1,1);
	low_set_seg_base(3,in2,1);
	low_int_mult(4,2,1,3,1);
	low_get_seg(4,out,LEN);
}

void random_op(unsigned int len)
{
	int i;
	for (i=0; i<len; i++)
	{
		out[i] = ecdh_get_random();
	}
}

void reg_copy()
{
	unsigned int i;
	for (i=0; i<aux; i++){
		out[i] = in1[i];
	}
}

void reg_xor()
{
	unsigned int i;
	for (i=0; i<aux; i++){
		out[i] = in1[i]^in2[i];
	}
}


void copy_acc_reg_3()
{
	unsigned int i;
	for (i=0; i<OPSIZE(LEN); i++){
		out[i] = ACCX[i];
		out[i+OPSIZE(LEN)] = ACCY[i];
		out[i+2*OPSIZE(LEN)] = ACCZ[i];
	}
}

void mod_import ()
{
	unsigned int L, i;

	L = in1[0];

	out[0] = 0;
	out[1] = 0;

	for(i = 1; i < L + 1; i++) out[i + 1] = in1[i];	
}

void standard_repr()
{
	low_set_seg_base(2,in1,LEN);
	low_standard_repr (4,2);
	low_get_seg(4,out,LEN);
}

void montgomery_repr()
{
	low_set_seg_base(2,in1,LEN);
	low_montgomery_repr(4,2,3);
	low_get_seg(4,out,LEN); 
}

void mod_red()
{
	low_set_seg_base(2,in1,LEN);
	low_mod_red(4,2);	
	low_get_seg(4,out,LEN);
}

void coord_invert()
{
	unsigned char winInv, digit, carry, j;
	short int i;
	unsigned int tmp,k;
	unsigned int *inInv, *acc1, *acc2, *exp;
	
	inInv = in1;
	
	//--- Define the window size ---

	// tmp is set to optimal window size according to LEN
	
	// window size is chosen such that the number 2^(w-1)+2 of required 
	// registers is at most the value in aux

	winInv = 3;

	if (LEN > 2)
	{
		if (aux >= (unsigned int)11) winInv = 4;
	}
	if (LEN > 7)
	{
		if (aux >= (unsigned int)19) winInv = 5;
	}

	//--- Define memory mapping ---

	// accumulator
	// two registers are required because of out-of-place multiplication
	acc1 = inInv + (1 << (winInv-1)) * OPSIZE(LEN); 
	acc2 = acc1 + OPSIZE(LEN);
	
	// exponent
	exp = acc2 + OPSIZE(LEN);
	
	//--- Precompute values x^3, x^5, x^7, etc. ---

	// Store x^2 at exp (temporary)
	in1 = inInv;
	out = exp;
	montgomery_square();

	// Compute x^(2i+1) = x^(2(i-1)+1).x^2
	for(i = 1; i < (1 << (winInv - 1)); i++)
	{
		in1 = inInv + (i-1)*OPSIZE(LEN);
		in2 = exp;
		out = inInv + i*OPSIZE(LEN);
		montgomery_mult();
	};

	//--- Compute exponent (modulus - 2) ---

	tmp = MOD[2];

	if(tmp == 1)
	{
		tmp = 0xFFFFFFFF;
		carry = 1;
	}
	else
	{
		tmp = tmp - 2;
		carry = 0;
	}
	
	exp[0] = tmp;
	for(k = 1; k < LEN; k++)
	{
		tmp = MOD[2+k];
		//tmp = tmp - carry;
		if(carry == 1)
		{
			if(tmp == 0)
			{
				tmp = 0xFFFFFFFF;
				carry = 1;
			}
			else {
				tmp = tmp - 1;
				carry = 0;
			}
		}
		exp[k] = tmp;
	}
	
	//--- Exponentiation loop ---
	
	// Copy input in acc1
	
	in1 = inInv;
	out = acc1;
	aux = OPSIZE(LEN);
	reg_copy();
	
	// Set i to the exponent msb -1
	
	// NB: non generic
	
	tmp = MOD[LEN + 1];
	i = 0;
	
	while(tmp > 0)
	{
		tmp = (unsigned int)tmp >> 1;
		i++;
	}
	
	i = 32 * (LEN - 1) + i - 2;
	
	while(i >= 0)
	{
		in1 = exp;
		in2 = (unsigned int*) ((unsigned int)i);
		tmp = get_bit();
		
		if(tmp == 0)
		{
			// Current exponent bit is 0
			// Square accumulator
			in1 = acc1;
			in2 = acc1;
			out = acc2;
			montgomery_mult();
			
			tmp = ((unsigned int)acc1);
			acc1 = acc2;
			acc2 = ((unsigned int*) tmp);
			
			// Decrement loop counter
			i--;
		}
		else {
			// Current exponent bit is 1
			// Get current digit value (in d) and bit-length (in aux)
			in1 = exp;
			in2 = (unsigned int*) ((unsigned int)i);
			aux = winInv;

			digit = get_digit_l2r();
			

			// Square accumulator aux times (aux is the digit bit-length)
			for(j = 0; j < aux; j++)
			{
				in1 = acc1;
				out = acc2;
				low_set_seg_base(2,in1,LEN);
				low_montgomery_mult(4,2,2);
				low_get_seg(4,out,LEN);

				tmp = ((unsigned int)acc1);
				acc1 = acc2;
				acc2 = ((unsigned int*) tmp);
			}
			
			// Multipliy accumulator with precomputed value x^d
			in1 = acc1;
			in2 = inInv + ((digit - 1)/2)*OPSIZE(LEN); 
			out = acc2;
			montgomery_mult();
			
			tmp = ((unsigned int)acc1);
			acc1 = acc2;
			acc2 = ((unsigned int*) tmp);
			
			// Decrement loop counter aux times
			i = i - aux;
		}
		
	}
	
	// copy result
	
	in1 = acc1;
	out = inInv;
	aux = OPSIZE(LEN);
	reg_copy();
	
}

unsigned char get_bit()
{
	return ((unsigned)in1[((unsigned)in2) >> 5] >> (((unsigned int)in2) & 0x1F)) & 1;
}

unsigned char get_digit_l2r()
{
	unsigned char indWord, indBit;
	unsigned int t0, t1;
	
	if(((unsigned int)in2) < aux)
	{
		t0 = in1[0] & ((1 << ((unsigned int)in2 + 1)) - 1);
		aux = ((unsigned int)in2) + 1;
	}
	
	else 
	{
		indWord = (unsigned long int)(((unsigned int)in2) + 1 - aux) >> 5;
		indBit = (((unsigned int)in2) + 1 - aux) & 0x1F;
		
		t0 = in1[indWord];
		t0 = (unsigned long int)t0 >> indBit;
		t0 = t0 & ((1 << aux) - 1);
		
		if(indBit + aux > 32)
		{
			t1 = in1[indWord + 1];
			t1 = t1 & ((1 << (((unsigned int)in2 + 1) & 0x1F)) - 1);
			t1 = t1 * (1 << (aux - (((unsigned int)in2 + 1) & 0x1F)));
			t0 = t0 + t1;
		}		
		
	}
	
	while((t0 & 0x1) == 0)
	{
		t0 = t0>>1; 
		aux--;
	}
	
	return t0;
}

unsigned int op_import()
{
  unsigned int L, i;

  L = in1[0];

  if (L > LEN + 1) L = LEN + 1;

  for(i = 0; i < L; i++)          	out[i] = in1[i + 1];
  for(i = L; i < OPSIZE(LEN); i++)  out[i] = 0;

  return L;
}

void op_export()
{
	aux = LEN;
	reg_copy();
}

void set_mod()
{
	low_set_mod(in1+2,aux);
}
