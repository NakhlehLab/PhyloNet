/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.containers;

import java.util.AbstractList;
import java.util.ArrayList;

public class ListOfArrays<Type> extends AbstractList<Type> implements
        IListOfArrays<Type> {

    private IArrayFactory<Type> arrayFactory;

    private ArrayList<IArray<Type>> arrays;

    /**
     * The number of elements currently stored in the arrays.
     */
    private int size;

    /**
     * The number of elements it is possible to store in the arrays if all
     * spaces are filled.
     */
    private int capacity;

    /**
     * The size of each array to be added.
     */
    private int growthIncrement;

    private static final int defaultGrowthIncrement = 100;
    private static final int defaultInitialCapacity = 0;

    public ListOfArrays(IArrayFactory<Type> arrayFactory) {
        init(arrayFactory, defaultGrowthIncrement, defaultInitialCapacity);
    }

    public ListOfArrays(IArrayFactory<Type> arrayFactory, int growthIncrement) {
        init(arrayFactory, growthIncrement, defaultInitialCapacity);
    }

    public ListOfArrays(IArrayFactory<Type> arrayFactory, int growthIncrement,
            int initialCapacity) {
        init(arrayFactory, growthIncrement, initialCapacity);
    }

    private void init(IArrayFactory<Type> arrayFactory, int growthIncrement,
            int initialCapacity) {

        // kliu - added IArray<Type> argument to fix compilation error
        arrays = new ArrayList<IArray<Type>>();

        assert (growthIncrement > 0);
        assert (initialCapacity > 0);

        this.arrayFactory = arrayFactory;
        this.size = 0;
        this.capacity = 0;
        this.growthIncrement = growthIncrement;

        while (capacity < initialCapacity) {
            enlarge();
        }
    }

    private void enlarge() {
        IArray<Type> next = arrayFactory.make(growthIncrement);
        arrays.add(next);

        capacity += growthIncrement;

        if (growthIncrement >= 10000) {
            System.out.println("Enlarging... capacity = " + capacity );
        }
    }

    private void shrink() {
        arrays.remove(arrays.size() - 1);

        capacity -= growthIncrement;
    }

    @Override
    public void add(int index, Type element) {
        // Check the bounds
        if (index < 0 || index > size()) {
            throw new IndexOutOfBoundsException("The index " + index
                    + " is invalid for a ListOfArrays of length " + size());
        }

        // Check if we need to enlarge
        if (size() == capacity) {
            enlarge();
        }

        // Shift everything to the right
        for (int i = size(); i > index; i--) {
            set(i, get(i - 1));
        }

        set(index, element);

        size++;
    }

    @Override
    public Type get(int index) {
        int arrayNum = index / growthIncrement;
        int indexInArray = index % growthIncrement;

        return arrays.get(arrayNum).get(indexInArray);
    }

    @Override
    public Type set(int index, Type element) {
        int arrayNum = index / growthIncrement;
        int indexInArray = index % growthIncrement;

        arrays.get(arrayNum).set(indexInArray, element);
        return element;
    }

    @Override
    public Type remove(int index) {
        // Check the bounds
        if (index < 0 || index >= size()) {
            throw new IndexOutOfBoundsException("The removal index " + index
                    + " is invalid for a ListOfArrays of length " + size()
                    + " (max index: " + (size() - 1) + ").");
        }

        Type result = get(index);

        // Shift everything to the left
        for (int i = index; i < size() - 1; i++) {
            set(i, get(i + 1));
        }

        // See if we can shrink
        if (size() < capacity - growthIncrement) {
            shrink();
        }

        return result;
    }

    @Override
    public int size() {
        return size;
    }
}
