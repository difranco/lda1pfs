package net.anthonydifranco.LDA;

import java.io.Serializable;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.DataOutputStream;

/**
 * Created by IntelliJ IDEA.
 * User: adifranco
 * Date: Jul 1, 2008
 * Time: 11:11:31 AM
 * To change this template use File | Settings | File Templates.
 */
public interface RandomNumberGenerator extends Serializable, Cloneable {
    boolean stateEquals(Object o);

    void readState(DataInputStream stream) throws IOException;

    void writeState(DataOutputStream stream) throws IOException;

    void setSeed(long seed);

    void setSeed(int[] array);

    boolean nextBoolean();

    boolean nextBoolean (float probability);

    boolean nextBoolean (double probability);

    int nextInt(int n);

    long nextLong(long n);

    double nextDouble();

    float nextFloat();

    void nextBytes(byte[] bytes);

    char nextChar();

    short nextShort();

    byte nextByte();

    double nextGaussian();
}
