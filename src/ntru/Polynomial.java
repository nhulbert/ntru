package ntru;

/**
 * Describes a truncated Polynomial with integer coefficients, and defines arithmetic operations of
 * these Polynomials over arbitrary modular rings.
 * 
 * @author Neil Hulbert
 */
public class Polynomial {
  
  /**
   * The coefficients of the polynomial, with the nth index pointing to the coefficient of the
   * nth power of the variable.
   */
  private int[] myData;

  /**
   * A constructor that takes a single argument specifying how many coefficients the polynomial has.
   * 
   * @param theSize the number of coefficients the polynomial has
   */
  Polynomial(final int theSize) {
    myData = new int[theSize];
  }
    
  /**
   * A constructor that encodes the given byte array into a polynomial of a given size.
   * 
   * @param theBytes the bytes to be encoded
   * @param theN the number of coefficients in the polynomial, must be greater-than-or
   *     equal-to theBytes*8
   */
  Polynomial(final byte[] theBytes, final int theN) {
    myData = new int[theN];
      
    int lim = theBytes.length * Byte.SIZE;
      
    for (int i = 0; i < lim && i < theN; i++) {
      myData[i] = (theBytes[i / Byte.SIZE] >>> (i % Byte.SIZE)) & 1;
    }
  }
  
  /**
   * A constructor that encodes the given integer array into a polynomial of the same length.
   * 
   * @param intArr the int array to be encoded
   */
  public Polynomial(int[] intArr) {
    myData = intArr.clone();
  }
    
  /**
   * A constructor that creates a new Polynomial equal to a previously existing Polynomial.
   * 
   * @param theOther the previously existing Polynomial to be copied from
   */
  public Polynomial(Polynomial theOther) {
    this(theOther.myData);
  }
    
  /**
   * Decodes a Polynomial into a byte array with bits set if the corresponding Polynomial
   * term is non-zero.
   * 
   * @return
   */
  public byte[] getBytes() {
    final byte[] bytes = new byte[(myData.length + Byte.SIZE - 1) / Byte.SIZE];
      
    for (int i = 0; i < myData.length; i++) {
      int value = 1;
      if (myData[i] == 0) {
        value = 0;
      }
          
      bytes[i / Byte.SIZE] ^= value << (i % Byte.SIZE);
    }
      
    return bytes;
  }
    
  /**
   * Reduces the Polynomial's coefficients in place to the smallest non-negative integers congruent
   * in the modular ring of integers with the given order.
   * 
   * @param mod the modulus to which the coefficients are to be reduced
   */
  public void reduceCoefsToMod(int mod) {
    for (int i = 0; i < myData.length; i++) {
      myData[i] = (myData[i] % mod + mod) % mod;
    }
  }
    
  /**
   * Reduces the Polynomial's coefficients in place to integers congruent in the modular ring of
   * integers with the given order, such that floor((1-mod)/2) <= n <= floor(mod/2) for any reduced
   * coefficient n.
   * 
   * @param mod the modulus to which the coefficients are to be reduced
   */
  public void convertToNegatives(int mod) {
    for (int i = 0; i < myData.length; i++) {
      myData[i] %= mod;
      if (myData[i] > mod / 2 || myData[i] < (1 - mod) / 2) {
        myData[i] -= (Math.abs(myData[i]) / myData[i]) * mod;
      }
    }
  }
    
  /**
   * Multiplies each coefficient of the Polynomial by a given integer value.
   * 
   * @param val the value by which each coefficient is to be multiplied
   * @param mod the modulus in which each coefficient's multiplication is to occur
   */
  public void multiplyScalarInPlace(int val, int mod) {
    for (int i = 0; i < myData.length; i++) {
      myData[i] *= val;
      myData[i] %= mod;
    }
  }

  /**
   * Adds to each coefficient of the Polynomial a given integer value.
   * 
   * @param theScalar the value by which each coefficient is to be incremented
   * @param theMod the modulus in which each coefficient's addition is to occur
   */
  public void addInPlace(int theScalar, int theMod) {
    for (int i = 0; i < myData.length; i++) {
      myData[i] += theScalar;
      myData[i] %= theMod;
    }
  }
    
  /**
   * Adds another Polynomial to the current one per coefficient and in place.
   * 
   * @param oth the Polynomial with which the current one is to be added
   * @param mod the modulus in which each coefficient's addition is to occur
   */
  public void addInPlace(Polynomial oth, int mod) {
    for (int i = 0; i < myData.length; i++) {
      myData[i] += oth.myData[i];
      myData[i] %= mod;
    }
  }

  /**
   * Subtracts another Polynomial from the current one per coefficient and in place.
   * 
   * @param oth the Polynomial to be subtracted from the current one
   * @param mod the modulus in which each coefficient's subtraction is to occur
   */
  public void subtractInPlace(Polynomial oth, int mod) {
    for (int i = 0; i < myData.length; i++) {
      myData[i] += -oth.myData[i] + mod;
      myData[i] %= mod;
    }
  }
    
  /**
   * Performs a convolution multiplication with another truncated polynomial of equal length,
   * returning a new Polynomial containing the result.
   * 
   * @param oth the Polynomial with which the current Polynomial is to be multiplied, must be
   *     the same length as the current Polynomial
   * @param mod the modulus in which each coefficient's multiplication/addition is to occur
   * @return a new Polynomial containing the result of the multiplication
   */
  public Polynomial getTruncatedMultiply(Polynomial oth, int mod) {
    Polynomial res = new Polynomial(myData.length);
      
    // There may well exist a better multiplication algorithm for sparse Polynomials:
    // this is intended to give approx. constant time given a fixed Polynomial length
    for (int i = 0; i < myData.length; i++) {
      res.myData[i] = myData[i] * oth.myData[0];
      for (int j = 1; j < myData.length; j++) {
        res.myData[i] += myData[(i + j) % myData.length] * oth.myData[myData.length - j];
        res.myData[i] %= mod;
      }
    }

    return res;
  }
    
  /**
   * Computes the multiplicative inverse of the Polynomial over the modular field of integers with
   * multiplicative inverses defined by the given array.
   * 
   * @param intInvTable a table holding the inverses of each non-unit, in order, within a
   *     modular multiplicative group, such that index n holds the multiplicative inverse
   *     of n+2 (mod l+1), with l being the length of intInvTable
   * @return a Polynomial that yields the Polynomial "1" when multiplied by the current Polynomial
   */
  public Polynomial inverse(int[] intInvTable) {
    int coefFieldOrder = intInvTable.length + 1;
        
    Polynomial[] polys = { new Polynomial(myData.length + 1), new Polynomial(myData) };
        
    int dividendInd = 0;
    int divisorInd = 1;
    
    Polynomial dividend = polys[dividendInd];
    Polynomial divisor = polys[divisorInd];
    
    // initialize the dividend to (x^N) - 1 which is
    // equivalent to 0 in the polynomial ring
    dividend.myData[0] = coefFieldOrder - 1;
    dividend.myData[myData.length] = 1;
        
    divisor.addInPlace(coefFieldOrder, coefFieldOrder);
        
    int[] coefInd = { myData.length, myData.length - 1 };

    // the Bezout coefficients, initialized to {0,1}
    Polynomial[] curInverse = { new Polynomial(myData.length), new Polynomial(myData.length) };
    curInverse[1].myData[0] = 1;        

    if (divisor.findNextNonzeroCoef(coefInd, divisorInd)) {
      return null;
    }
        
    // The loop for the Euclidean algorithm for truncated polynomials
    while (coefInd[divisorInd] > 0) {
      Polynomial quotient = new Polynomial(myData.length);
            
      // The loop for Polynomial long division
      while (coefInd[dividendInd] >= coefInd[divisorInd]) {
        if (dividend.myData[coefInd[dividendInd]] != 0) {
          int quotientPower = coefInd[dividendInd] - coefInd[divisorInd];
          int quotientCoef = (dividend.myData[coefInd[dividendInd]]
              * intInvTable[divisor.myData[coefInd[divisorInd]] - 1])
              % coefFieldOrder;
          quotient.myData[quotientPower] = quotientCoef;

          // the definition of quotientCoef ensures that the leading dividend coefficient is always
          // eliminated
          dividend.myData[coefInd[dividendInd]] = 0; 
          for (int i = coefInd[dividendInd] - 1; i >= quotientPower; i--) {
            dividend.myData[i] += coefFieldOrder
                - ((quotientCoef * divisor.myData[i - quotientPower])
                    % coefFieldOrder);
            dividend.myData[i] %= coefFieldOrder;
          }
        }
        
        coefInd[dividendInd]--;
      }
      
      if (dividend.findNextNonzeroCoef(coefInd, dividendInd)) {
        return null;
      }

      Polynomial update = quotient.getTruncatedMultiply(curInverse[divisorInd], coefFieldOrder);
      update.multiplyScalarInPlace(coefFieldOrder - 1, coefFieldOrder); // multiply by -1
      curInverse[dividendInd].addInPlace(update, coefFieldOrder);
            
      divisorInd = dividendInd;
      dividendInd = 1 - dividendInd;

      dividend = polys[dividendInd];
      divisor = polys[divisorInd];
    }

    curInverse[divisorInd].multiplyScalarInPlace(intInvTable[divisor.myData[0] - 1],
        coefFieldOrder);

    return curInverse[divisorInd];
  }

  /**
   * Finds the biggest non-zero coefficient at an index no greater than coefInd[polyInd] within the
   * current Polynomial. Intended so the current Polynomial can be prepared to be used as a divisor
   * in Polynomial long division.
   * 
   * @param coefInd an array containing indices of Polynomials
   * @param polyInd the index of coefInd corresponding to the Polynomial
   * @return a boolean value representing whether such a non-zero coefficient exists in the current
   *     Polynomial.
   */
  private boolean findNextNonzeroCoef(int[] coefInd, int polyInd) {
    while (myData[coefInd[polyInd]] == 0) {
      coefInd[polyInd]--;
      if (coefInd[polyInd] < 0) {
        return true;
        // if the last non-zero remainder's degree is > 0, there is no
        // inverse and the Euclidean algorithm will terminate
      }
    }
      
    return false;
  }
    
  /**
   * Returns the length of the current Polynomial, or equivalently, what power of the variable
   * is equivalent to 1 in the truncated Polynomial ring.
   * 
   * @return the length of the current Polynomial
   */
  public int length() {
    return myData.length;
  }
    
  /**
   * Returns a full, human-readable representation of the Polynomial as a String.
   */
  public String toString() {
    StringBuilder sb = new StringBuilder();

    for (int i = 0; i < myData.length; i++) {
      if (myData[i] != 0) {
        if (myData[i] > 0) {
          sb.append("+");
        } else {
          sb.append("-");
        }
        if (i == 0 || Math.abs(myData[i]) > 1) {
          sb.append(Math.abs(myData[i]));
        }
        if (i > 0) {
          sb.append("x");
          if (i > 1) {
            sb.append("^" + i);
          }
        }
      }
    }

    if (sb.length() == 0) {
      return "0";
    }
    if (sb.charAt(0) == '+') {
      return sb.substring(1);
    }
    return sb.toString();
  }
    
  /**
   * Returns whether this Polynomial is equal to another Object, up to Polynomial length and
   * per-coefficient integer equality (not modular equivalence).
   */
  @Override
  public boolean equals(Object oth) {
    if (this == oth) {
      return true;
    }
      
    if (oth == null) {
      return false;
    }
      
    if (oth.getClass() != getClass()) {
      return false;
    }
      
    Polynomial othPoly = (Polynomial)oth;
      
    if (myData.length != othPoly.myData.length) {
      return false;
    }
      
    for (int i = 0; i < myData.length; i++) {
      if (myData[i] != othPoly.myData[i]) {
        return false;
      }
    }
      
    return true;
  }
}