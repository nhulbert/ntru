package ntru;

import java.security.SecureRandom;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

/**
 * A class defining an NTRUEncrypt key pair, which includes the private key comprised of Polynomials
 * f and fpInv, as well as the public key comprised of the Polynomial h. Also contains fields
 * representing the parameters N, P, Q, df, dg, and dr of the key system in which the key pair
 * operates.
 * 
 * <p>Based on the paper NTRU: A Ring-Based Public Key Cryptosystem
 * by Jeffrey Hoffstein, Jill Pipher, and Joseph H. Silverman
 * 
 * @author Neil Hulbert
 */
public class KeyPair {
  /**
   * half of the private key, a Polynomial that must be invertible mod modQ and modP.
   */
  private Polynomial polyF;
    
  /**
   * half of the private key, a Polynomial that is the inverse of polyF mod modP.
   */
  private Polynomial polyFpInv;
  
  /**
   * the Polynomial that represents the public key, must equal the inverse of polyF
   * mod lengthN multiplied successively by modP and then by some "small" Polynomial g.
   */
  private Polynomial polyH;
  
  /**
   * the number of coefficients in a Polynomial in the key system, equivalently,
   * the power of the variable which is equivalent to 1 in the truncated Polynomial ring.
   */
  private int lengthN;
  
  /**
   * the large integer modulus of the key system, must be relatively prime with modP.
   */
  private int modQ;
  
  /**
   * the small integer modulus of the key system, must be relatively prime with q.
   */
  private int modP;
  
  /**
   * the number of coefficients equal to each 1 and one greater than the number of
   * coefficients equal to -1 within Polynomial polyF.
   */
  private int df;
  
  /**
   * the number of coefficients equal to each 1 and -1 in the Polynomial g.
   */
  private int dg;
  
  /**
   * the number of coefficients equal to each 1 and -1 in the Polynomial r, the
   * Polynomial which acts as the blinding value for the message during encryption.
   */
  private int dr;
  
  /**
   * The constructor for creating an NTRU KeyPair.
   * 
   * @param f half of the private key, a Polynomial that must be invertible mod n and p
   * @param fpInv half of the private key, a Polynomial that is the inverse of f mod p
   * @param h the Polynomial that represents the public key, must equal the inverse of f
   *     mod n multiplied successively by p and then by some "small" Polynomial g
   * @param n the number of coefficients in a Polynomial in the key system, equivalently,
   *     the power of the variable which is equivalent to 1 in the truncated Polynomial ring
   * @param p the small integer modulus of the key system, must be relatively prime with q
   * @param q the large integer modulus of the key system, must be relatively prime with p
   * @param df the number of coefficients equal to each 1 and one greater than the number of
   *     coefficients equal to -1 within Polynomial f
   * @param dg the number of coefficients equal to each 1 and -1 in the Polynomial g
   * @param dr the number of coefficients equal to each 1 and -1 in the Polynomial r, the
   *     Polynomial which acts as the blinding value for the message during encryption
   */
  public KeyPair(Polynomial f, Polynomial fpInv, Polynomial h, int n,
                int p, int q, int df, int dg, int dr) {
    this.polyF = f;
    this.polyFpInv = fpInv;
    this.polyH = h;
    this.lengthN = n;
    this.modP = p;
    this.modQ = q;
    this.df = df;
    this.dg = dg;
    this.dr = dr;
  }
  
  /**
   * Decrypts with this KeyPair the passed in Polynomial which represents an NTRU encrypted message.
   * 
   * @param encryptedMessage the Polynomial representing the NTRU encrypted message
   * @return the decrypted message as a byte array
   */
  public byte[] decrypt(Polynomial encryptedMessage) {
    Polynomial a = polyF.getTruncatedMultiply(encryptedMessage, modQ);
    a.convertToNegatives(modQ);
    a.convertToNegatives(modP);
    
    Polynomial decryptedMessage = a.getTruncatedMultiply(polyFpInv, modP);
    decryptedMessage.convertToNegatives(modP);
    
    return decryptedMessage.getBytes();
  }
  
  /**
   * Returns the Polynomial representing the public key of the key pair.
   * 
   * @return the Polynomial representing the public key of the key pair
   */
  public Polynomial getPublicKey() {
    return polyH;
  }
  
  /**
   * Encrypts a byte array representing the secret message.
   * 
   * @param theMessage the byte array representing the secret message
   * @param h the Polynomial representing the public key
   * @param n the integer representing the length of Polynonials within the key system,
   *     alternatively, the power of the variable equivalent to 1 in the truncated Polynomial ring
   * @param q the large integer modulus of the key system
   * @param dr the number of coefficients equal to each 1 and -1 to be included in the blinding
   *     Polynomial r
   * @param random the random number generator used to generate the blinding Polynomial r
   * @return a Polynomial that represents the encrypted message
   */
  public static Polynomial encrypt(byte[] theMessage, Polynomial h, int n,
                                   int q, int dr, Random random) {
    Polynomial messagePoly = new Polynomial(theMessage, n);
    
    Polynomial obscurer = generateRandomPolynomial(n, dr, dr, random);
    Polynomial encryptedMessage = h.getTruncatedMultiply(obscurer, q);
    encryptedMessage.addInPlace(messagePoly, q);
    
    return encryptedMessage;
  }
  
  /**
   * Randomly generates a key pair within the key system specified by the parameters.
   * 
   * @param theN the integer representing the length of Polynomials within the key system,
   *     alternatively, the power of the variable equivalent to 1 in the truncated Polynomial ring
   * @param theP the small integer modulus of the key system
   * @param theQ the large integer modulus of the key system
   * @param df the number of coefficients equal to each 1 and one greater than the number of
   *     coefficients equal to -1 to be included in the generated Polynomial f
   * @param dg the number of coefficients equal to each 1 and -1 to be included in the generated
   *     Polynomial g
   * @param dr the number of coefficients equal to each 1 and -1 to be included in the generated
   *     blinding Polynomial r
   * @param random the source of randomness to be used in generating the key pair
   * @return the key pair holding both the private and public keys
   */
  public static KeyPair generateKeyPair(int theN, int theP, int theQ, int df,
                                        int dg, int dr, Random random) {
    Map<Integer, int[]> intInversesMap = new HashMap<>();
      
    Map<Integer, Integer> factorsOfQ = intFactor(theQ);
    Map<Integer, Integer> factorsOfP = intFactor(theP);
      
    Stream.of(factorsOfQ, factorsOfP).forEach(factors -> 
        factors.forEach((factor, freq) ->
            intInversesMap.putIfAbsent(factor, generateInverses(factor))));
        
    Polynomial g = generateRandomPolynomial(theN, dg, dg, random);
    Polynomial f = null;
    Polynomial fqInv = null;
    Polynomial fpInv = null;
        
    while (fpInv == null || fqInv == null) {
      // number of 1's and -1's must differ in f, otherwise it won't be invertible
      f = generateRandomPolynomial(theN, df, df - 1, random);
      fpInv = compositeModInverse(f, factorsOfP, intInversesMap);
      fqInv = compositeModInverse(f, factorsOfQ, intInversesMap);
    }
        
    fqInv.multiplyScalarInPlace(theP, theQ);
    Polynomial h = fqInv.getTruncatedMultiply(g, theQ);
        
    return new KeyPair(f, fpInv, h, theN, theP, theQ, df, dg, dr);
  }
    
  /**
   * Generates a random Polynomial from the set of Polynomials specified in the parameters.
   * 
   * @param length the length of the Polynomial, alternatively, the power of the variable
   *     that is equivalent to 1 in the truncated Polynomial ring
   * @param numPosOnes the number of 1's that should appear as coefficients in the generated
   *     Polynomial, must be less than or equal to length - numNegOnes
   * @param numNegOnes the number of -1's that should appear as coefficients in the generated
   *     Polynomial, must be less than or equal to length - numPosOnes
   * @param random the source of randomness from which the Polynomial should be generated
   * @return a randomly generated Polynomial
   */
  private static Polynomial generateRandomPolynomial(int length, int numPosOnes,
                                                       int numNegOnes, Random random) {
    int[] data = new int[length];
        
    int posOnesCount = 0;
    int negOnesCount = 0;
        
    for (int totalEvaluated = 0; totalEvaluated < data.length; totalEvaluated++) {
      int testVal = random.nextInt(data.length - totalEvaluated);
      if (testVal < numPosOnes - posOnesCount) {
        data[totalEvaluated] = 1;
        posOnesCount++;
      } else if (testVal < numPosOnes + numNegOnes - posOnesCount - negOnesCount) {
        data[totalEvaluated] = -1;
        negOnesCount++;
      }
    }
        
    return new Polynomial(data);
  }
    
  /**
   * Generates an array which holds, in order, the multiplicative inverses of all units besides
   * 1 in a modular ring.
   * 
   * @param mod the order of the modular ring
   * @return an array which holds, in order, the multiplicative inverses of all units besides
   *     1 in Z_mod
   */
  private static int[] generateInverses(int mod) {
    int[] res = new int[mod - 1];
    res[0] = 1;
      
    for (int i = 2; i < mod; i++) {
      if (res[i - 1] == 0) {
        res[i - 1] = intInverse(i, mod);
        if (res[i - 1] != 1) {
          res[res[i - 1] - 1] = i;
        }
      }
    }
      
    return res;
  }
     
  /**
   * Computes the inverse of an element; returns 1 if an inverse does not exist.
   * 
   * @param val the element to be inverted
   * @param mod the modulus in which the inverse is to be found
   * @return the inverse of val in modulus mod, or 1 if no such inverse exists
   */
  private static int intInverse(int val, int mod) {
    return (intExtendedEuclidean(val, mod)[0] % mod + mod) % mod;
  }
     
  /**
   * Performs the extended Euclidean algorithm to find what linear combination of 
   * the two integer arguments yields 1.
   * 
   * @param a the first integer value
   * @param b the second integer value
   * @return a two-element array of integers p such that p[0]*a + p[1]*b = 1
   */
  private static int[] intExtendedEuclidean(int a, int b) {
    int[] vals = {a, b};
    int[][] bezout = {{1, 0}, {0, 1}};
   
    int ind = (a > b) ? 1 : 0;
   
    while (vals[ind] != 0) {
      int quotient = vals[1 - ind] / vals[ind];
      for (int i = 0; i < 2; i++) {
        bezout[1 - ind][i] -= quotient * bezout[ind][i];
      }
      vals[1 - ind] %= vals[ind];
     
      ind = 1 - ind;
    }
   
    return bezout[1 - ind];
  }
     
  /**
   * Generates a Map representing the prime factorization of a given positive integer value.
   * 
   * @param p a positive integer to be factored into its primes
   * @return a Map holding integer keys representing each prime factor of p, with each
   *     respective value corresponding to the multiplicity of the prime in p's prime factorization
   */
  private static Map<Integer, Integer> intFactor(int p) {       
    Map<Integer, Integer> result = new HashMap<>();
    
    for (int i = 2; i * i <= p; i++) {
      while (p % i == 0) {
        result.put(i, result.getOrDefault(i, 0) + 1);
        p /= i;
      }
    }
    
    if (p != 1) {
      result.put(p, result.getOrDefault(p, 0) + 1);
    }
    
    return result;
  }
  
  /**
   * Computes the inverse of a truncated Polynomial over a modular field having a possibly composite
   * order. First computes the inverses mod each prime factor j of the desired modulus, then lifts
   * each inverse to mod j^k (k being some positive integer), and finally combines these inverses to
   * achieve the final integer modulus for the coefficients.
   * 
   * @param p the truncated Polynomial to be inverted
   * @param factors a Map representing the prime factorization of the final desired modulus over
   *     which the inverse is to be taken
   * @param intInversesMap a Map holding all multiplicative inverses in the various modular integer
   *     rings which will be needed to compute the Polynomial's inverse
   * 
   * @return a truncated Polynomial that when multiplied by p, yields 1 in the integer modulus
   *     described by the factors parameter
   */
  private static Polynomial compositeModInverse(Polynomial p, Map<Integer, Integer> factors,
                                                Map<Integer, int[]> intInversesMap) {       
    Polynomial result = factors.entrySet().parallelStream()
        .map(e -> liftInverse(p, e.getKey(), e.getValue(), intInversesMap))
        .reduce(new PolynomialMod(p.length(), 1), PolynomialMod::combine).poly;
       
    return result;
  }
     
  /**
   * Finds the inverse of a truncated Polynomial, in a given integer modulus that is a
   * positive integer power of a prime. This uses a formula derived from Hensel's lemma.
   * 
   * @param p the Polynomial to be inverted
   * @param mod the base of the exponential expression describing the final desired modulus, must
   *     be prime
   * @param power the power of the exponential expression describing the final desired modulus
   * @param intInversesMap a Map holding all multiplicative inverses in the various modular
   *     integer rings which will be needed to compute the Polynomial's inverse
   * @return a truncated Polynomial that when multiplied by p, yields 1 in the integer modulus
   *     mod^power
   */
  private static PolynomialMod liftInverse(Polynomial p, int mod, int power,
                                           Map<Integer, int[]> intInversesMap) {
    Polynomial current = p.inverse(intInversesMap.get(mod));
       
    int curPower = 1;
    int curMod = mod;
       
    while (curPower < power) {
      curMod *= mod;
      curPower++;
         
      if (current != null) {
        Polynomial temp = current.getTruncatedMultiply(current, curMod)
            .getTruncatedMultiply(p, curMod);
        current.multiplyScalarInPlace(2, curMod);
        current.subtractInPlace(temp, curMod);         
      }
    }
       
    return new PolynomialMod(current, curMod);
  }
     
  /**
   * A class that describes an element of the set of Polynomials with integer coefficients
   * crossed with the integers, with the integer representing the order of the modular
   * integer ring over which the Polynomial is taken.
   * 
   * @author Neil Hulbert
   */
  private static class PolynomialMod {
    /**
     * The truncated polynomial.
     */
    Polynomial poly;
       
    /**
     * The integer modulus of the polynomial's coefficients.
     */
    int mod;
       
    /**
     * Creates a new PolynomialMod Object given the initial values of both fields.
     * 
     * @param poly the initial value of the Polynomial
     * @param mod the integer modulus of the polynomial's coefficients
     */
    public PolynomialMod(Polynomial poly, int mod) {
      this.poly = poly;
      this.mod = mod;
    }
       
    /**
     * Creates a new PolynomialMod object given the length of the Polynomial and the
     * initial integer modulus.
     * 
     * @param length the length of the Polynomial, to be initialized with all zeros
     * @param mod the initial integer modulus
     */
    public PolynomialMod(int length, int mod) {
      this(new Polynomial(length), mod);
    }
       
    /**
     * Combines two PolynomialMods with moduli a and b such that the resulting PolynomialMod's new
     * integer modulo is a*b and that its Polynomial is unchanged taken over the original integer
     * moduli. It achieves this through application of the Chinese Remainder Theorem applied to the
     * truncated Polynomial rings.
     * 
     * @param p1 the first PolynomialMod to be combined
     * @param p2 the second PolynomialMod to be combined
     * @return the PolynomialMod resulting from combining p1 and p2
     */
    public static PolynomialMod combine(PolynomialMod p1, PolynomialMod p2) {
      int newMod = p1.mod * p2.mod;
         
      if (p1.poly == null || p2.poly == null) {
        return new PolynomialMod(null, newMod);
      }
         
      int[] bezout = intExtendedEuclidean(p1.mod, p2.mod);
      
      bezout[0] = (bezout[0] + newMod) % newMod;
      bezout[1] = (bezout[1] + newMod) % newMod;
      
      Polynomial result = new Polynomial(p1.poly);
      result.multiplyScalarInPlace(bezout[1] * p2.mod, newMod);
       
      Polynomial temp = new Polynomial(p2.poly);
      temp.multiplyScalarInPlace(bezout[0] * p1.mod, newMod);
       
      result.addInPlace(temp, newMod);
       
      return new PolynomialMod(result, newMod);
    }
  }
     
  /**
   * Returns n, the number of coefficients in a Polynomial in the key system, equivalently,
   * the power of the variable which is equivalent to 1 in the truncated Polynomial ring.
   * 
   * @return an integer representing the number of coefficients in a Polynomial in the key
   *     system, or equivalently, the power of the variable which is equivalent to 1 in the
   *     truncated Polynomial ring
   */
  public int getN() {
    return lengthN;
  }

  /**
   * Returns q, the large integer modulus of the key system. 
   * 
   * @return an integer value representing the large integer modulus of the key system
   */
  public int getQ() {
    return modQ;
  }

  /**
   * Returns p, the small integer modulus of the key system.
   * 
   * @return an integer value representing the small integer modulus of the key system
   */
  public int getP() {
    return modP;
  }

  /**
   * Returns df, which is both the number of 1's in the Polynomial f, and one
   * greater than the number of -1's in the Polynomial f.
   * 
   * @return an integer value representing df, which is both the number of 1's in the Polynomial f,
   *     and one greater than the number of -1's in the Polynomial f
   */
  public int getDf() {
    return df;
  }

  /**
   * Returns dg, which is both the number of 1's in the Polynomial g, and the number of -1's in the
   * Polynomial g.
   * 
   * @return an integer value representing dg, which is both the number of 1's in the Polynomial g,
   *     and the number of -1's in the Polynomial g
   */
  public int getDg() {
    return dg;
  }

  /**
   * Returns dr, which is both the number of 1's in the blinding Polynomial r, and the number of
   * -1's in the Polynomial r.
   * 
   * @return an integer value representing dr, which is both the number of 1's in the blinding
   *     Polynomial r, and the number of -1's in the Polynomial r
   */
  public int getDr() {
    return dr;
  }

  /**
   * The entry point to a demo of the NTRUEncrypt public key system.
   * 
   * @param args the command line arguments, which are ignored
   */
  public static void main(String[] args) {
    Random random = new SecureRandom();
     
    KeyPair pair = KeyPair.generateKeyPair(503, 2, 253, 155, 100, 65, random);
    
    Polynomial encrypted1 = KeyPair.encrypt("Hello World!".getBytes(), pair.getPublicKey(),
                                           pair.getN(), pair.getQ(), pair.getDr(), random);
    Polynomial encrypted2 = KeyPair.encrypt("The quick brown fox, etc.".getBytes(),
                                            pair.getPublicKey(), pair.getN(), pair.getQ(),
                                            pair.getDr(), random);
     
    System.out.println(encrypted1 + "\n" + encrypted2 + "\n");
     
    // implementation does not yet include padding and
    // will not be secure in all applications
    String decrypted1 = new String(pair.decrypt(encrypted1)); 
    String decrypted2 = new String(pair.decrypt(encrypted2));
    
    // Strings are longer than original message: they are padded with the null character
    System.out.println(decrypted1); 
    System.out.println(decrypted2);
  }
}