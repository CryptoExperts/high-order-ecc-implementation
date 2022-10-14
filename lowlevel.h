// -------------------------------------------------------------------------- //
//---- Language:         C
//---- Module name:      lowlevel.h
//---- Developed:        Sonia Belaïd & Matthieu Rivain
//---- Contact:          {sonia.belaid,matthieu.rivain}@cryptoexperts.com
//---- Copyright:        CryptoExperts
// -------------------------------------------------------------------------- //

//--- TRNG ---//

void 			RNG_on(void);
void 			RNG_off(void);
unsigned int 	random_word(void);

//--- Low-level Functions ---//

/**
 * Switch on/off the coprocessor
 */
void copro_on(void);
void copro_off(void);

/**
 * Load the modulus in the coprocessor
 * - mod: modulus
 * - len: lenght of modulus
 */
void low_set_mod(unsigned int *mod, unsigned int len);

/**
 * Record in memory
 * - src: source register
 * - dst: index of the memory segment
 * - len: length to copy
 */
void low_set_seg_base(unsigned int dst, unsigned int *src, unsigned int len);
/**
 * Record in memory and add zeroes when shorter than LEN
 * - src: source register
 * - dst: index of the memory segment
 * - len: length to copy
 */
void low_set_seg(unsigned int dst, unsigned int *src, unsigned int len);
/**
 * Load from memory
 * - src: index of the memory segment
 * - dst: destination register
 * - len: length to copy
 */
void low_get_seg(unsigned int src, unsigned int *dst, unsigned int len);
/**
 * Copy in memory
 * - src: index of the source memory segment
 * - dst: index of the destination memory segment
 */
void low_copy_seg(unsigned int dst, unsigned int src);
/**
 * Compute the Montgomery inverse 
 * - tmp: index of a temporary register
 * - res: index of the Montgomery inverse
 */
void low_montgomery_inv(unsigned int res, unsigned int tmp);
/**
 * Compute the Montgomery representation
 * - tmp: index of a temporary register
 * - op: index of the source register
 * - res: index of the destination register
 */
void low_montgomery_repr(unsigned int res, unsigned int op, unsigned int tmp);
/**
 * Recover the standard representation
 * - op: index of the source register
 * - res: index of the destination register
 */
void low_standard_repr(unsigned int res, unsigned int op);
/**
 * Compute the Montgomery multiplication between op1 and op2
 * without setting the base register
 * - op1: index of the first source register
 * - op2: index of the second source register
 * - res: index of the destination register
 */
void low_montgomery_mult_after(unsigned int res, unsigned int op1, unsigned int op2);
/**
 * Compute the Montgomery multiplication between op1 and op2
 * - op1: index of the first source register
 * - op2: index of the second source register
 * - res: index of the destination register
 */
void low_montgomery_mult(unsigned int res, unsigned int op1, unsigned int op2);
/**
 * Compute the multiplication between op1 and op2
 * - op1: index of the first source register
 * - len1: length of the first input
 * - op2: index of the second source register
 * - len2: length of the second input
 * - res: index of the destination register
 */
void low_int_mult(unsigned int res, unsigned int op1, unsigned int len1, unsigned int op2, unsigned int len2);
/**
 * Compute the arithmetic operation op1+(op2*op3)
 * - op1: index of the first source register
 * - len1: length of the first input
 * - op2: index of the second source register
 * - len2: length of the second input
 * - op3: index of the third source register
 * - len3: length of the third input
 * - res: index of the destination register
 */
void low_int_macc(unsigned int res, unsigned int op1, int len1, unsigned int op2, int len2, unsigned int op3, int len3);
/**
 * Compute the arithmetic subtraction modulo MOD (op1-op2)
 * - op1: index of the first source register
 * - op2: index of the second source register
 * - res: index of the destination register
 */
void low_mod_sub(unsigned int res, unsigned int op1, unsigned int op2);
/**
 * Compute the opposite (MOD - op1)
 * - op1: index of the source register
 * - res: index of the destination register
 */
void low_mod_sub_mod(unsigned int res, unsigned int op1);
/**
 * Compute the arithmetic addition modulo MOD
 * - op1: index of the first source register
 * - op2: index of the second source register
 * - res: index of the destination register
 */
void low_mod_add(unsigned int res, unsigned int op1, unsigned int op2);
/**
 * Compute the arithmetic doubling modulo MOD
 * - op1: index of the source register
 * - res: index of the destination register
 */
void low_mod_doubl(unsigned int res, unsigned int op1);
/**
 * Compute the modular reduction modulo MOD
 * - op1: index of the source register
 * - res: index of the destination register
 */
void low_mod_red(unsigned int res, unsigned int op1);




