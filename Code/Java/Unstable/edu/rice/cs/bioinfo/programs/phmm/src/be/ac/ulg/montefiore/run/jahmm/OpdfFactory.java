/* jahmm package - v0.6.1 */

/*
 *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package be.ac.ulg.montefiore.run.jahmm;


/**
 * Classes implementing this interface are able to generate observation
 * probability distribution functions.
 * The classes implementing <code>OpdfFactory</code> are supposed to generate
 * only a certain kind of distribution (e.g. Gaussian).
 */
public interface OpdfFactory<D extends Opdf<?>>
{    
    /**
     * Generates a new observation probability distribution function.
     * 
     * @return The new opdf.
     */
    public D factor();
}
