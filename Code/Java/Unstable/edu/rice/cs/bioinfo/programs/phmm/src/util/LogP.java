/**
 * log(x), where -\infty supported,
 * i.e. modified to support log(0)
 * all ops base e.
 *
 * set to 0.0 = log(1.0) initially.
 */

package util;

public class LogP {
    protected double d;
    protected boolean negativeInfinityFlag;

    protected static final String NEGATIVE_INFINITY_STRING = "NEGATIVE_INFINITY";
    public static final double EXPONENT_DIFFERENCE_ZERO_DELTA = 1e-30;

    // for efficiency - don't mutate these!
    public static final LogP ZERO_LOG_P = new LogP((double) 0.0);
    public static final LogP ONE_LOG_P = new LogP((double) 1.0);
    public static final LogP TWO_LOG_P = new LogP((double) 2.0);

    public LogP () {
	this(0.0);
    }

    public LogP (double x) {
	set(x);
    }

    public LogP (LogP m) {
	set(m);
    }

    /**
     * Directly set the members of this class.
     * Warning - member d is log of value!!!
     */
    public LogP (double inD, boolean inNegativeInfinityFlag) {
	set(inD, inNegativeInfinityFlag);
    }

    /**
     * Directly set the members of this class.
     * Warning - member d is log of value!!!
     */
    public void set (double inD, boolean inNegativeInfinityFlag) {
	this.d = inD;
	this.negativeInfinityFlag = inNegativeInfinityFlag;
    }

    public void set (double x) {
	if (x <= 0.0) {
	    d = 0.0;
	    negativeInfinityFlag = true;
	}
	else {
	    d = Math.log(x);
	    negativeInfinityFlag = false;
	}
    }

    public void set (LogP m) {
	if (m == null) {
	    System.err.println ("ERROR: cannot set LogP using a null reference! Returning.");
	    return;
	}

	this.d = m.d;
	this.negativeInfinityFlag = m.negativeInfinityFlag;
    }

    /**
     * WARNING: 
     * don't call this if you expect checkNegativeInfinity()
     * to return true.
     */
    public double getExponent () {
	if (checkNegativeInfinity()) {
	    return (Double.NEGATIVE_INFINITY);
	}
	else {
	    return (this.d);
	}
    }

    /**
     * WARNING: 
     * Beware of overflow/underflow!
     */
    public double doubleValue () {
	if (checkNegativeInfinity()) {
	    return (0.0);
	}
	else {
	    return (Math.exp(d));
	}
    }

    public boolean checkNegativeInfinity () {
	return (negativeInfinityFlag);
    }

    /**
     * Use a slightly faster formula.
     * rhat = ln (p + q) = phat + ln ( 1 + e^{qhat - phat} )
     * where 
     * qhat = ln q
     * phat = ln p
     *
     * Set subtractFlag to 1 to subtract instead.
     *
     * WARNING: the log operation below may return negative infinity
     * if the subtraction operation evaluates to a negative result.
     * No guard against that here in this code. Caller's
     * responsibility to guard against that.
     */
    public void add (LogP m, boolean subtractFlag) {
	if (m.checkNegativeInfinity()) {
	    return;
	}
	else if (this.checkNegativeInfinity()) {
	    this.set(m);
	}
	// need to handle subtract == 0 properly
	// better to guard by difference close to zero by super-small
	// amount
	// do we really care about super-small amount of precison anyways??
	// meh...
	// just lose the contribution - not important, most likely
	else if (subtractFlag && (Math.abs(m.d - this.d) < EXPONENT_DIFFERENCE_ZERO_DELTA)) {
	    this.set(0.0);
	}
	else {
	    if (subtractFlag) {
		this.d = this.d + Math.log(1.0 - Math.exp(m.d - this.d));
	    }
	    else {
		this.d = this.d + Math.log(1.0 + Math.exp(m.d - this.d));
	    }
	}
    }

    public void add (LogP m) {
	add (m, false);
    }

    public void subtract (LogP m) {
	add (m, true);
    }

    public void power (double k) {
	if (this.checkNegativeInfinity()) {
	    return;
	}
	else {
	    this.d *= k;
	}
    }

    /**
     * Set divideFlag to true for divide, otherwise multiply.
     */
    public void multiply (LogP m, boolean divideFlag) {
	if (m.checkNegativeInfinity()) {
	    this.d = 0;
	    this.negativeInfinityFlag = true;
	}
	else if (this.checkNegativeInfinity()) {
	    return;
	}
	else {
	    if (divideFlag) {
		this.d -= m.d;
	    }
	    else {
		this.d += m.d;
	    }
	}
    }

    public void multiply (LogP m) {
	multiply(m, false);
    }

    public void divide (LogP m) {
	multiply(m, true);
    }
    

    public boolean isGreaterThan(LogP m) {
	if (this.checkNegativeInfinity() && m.checkNegativeInfinity()) {
	    return (false);
	}
	else if (this.checkNegativeInfinity() && !m.checkNegativeInfinity()) {
	    return (false);
	}
	else if (!this.checkNegativeInfinity() && m.checkNegativeInfinity()) {
	    return (true);
	}
	else {
	    return (this.d > m.d);
	}
    }

    public String toString () {
	return (toString(false));
    }

    /**
     * Set expFlag to true to exponentiate logs.
     */
    public String toString (boolean expFlag) {
	if (checkNegativeInfinity()) {
	    if (expFlag) {
		return (Integer.toString(0));
	    }
	    else {
		return (NEGATIVE_INFINITY_STRING);
	    }
	}
	else {
	    if (expFlag) {
		// beware overflow!
		return (Double.toString(Math.exp(d)));
	    }
	    else {
		return (Double.toString(d));
	    }
	}
    }


}
