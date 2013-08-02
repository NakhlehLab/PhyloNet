/**
 * UnivariateFunction with exposed methods to work with a single mutable parameter.
 */

package optimize;

import org.apache.commons.math3.analysis.UnivariateFunction;

public interface MutableUnivariateFunction<P extends Parameter> extends UnivariateFunction {
    // e.g., branch length constructor will take in branch reference
    // subsequently deal with length get/set

    /**
     * Set the value of the mutable parameter.
     */
    public void setParameter (P p);

    /**
     * Get the value of the mutable parameter.
     */
    public P getParameter ();
}
