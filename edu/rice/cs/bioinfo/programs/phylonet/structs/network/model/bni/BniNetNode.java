/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;

import java.util.*;

/**
 * This class implements the methodes declared in the interface NetNode.
 *
 * @author Cuong Than
 *
 * @param <T> indicates the type of additional data this node stores.f
 *
 * June 12, 06: Created
 */
public class BniNetNode<T> implements NetNode<T> {
	/**
	 * This constructor instantiates an isolated, no-name node.
	 */
	public BniNetNode()
	{
		_data = null;
		_name = NO_NAME;
		_children = null;
		_parents = null;
		_parent_distances = null;
        _parent_probabiliites = null;
        _parent_support = null;
		_parent_gen_times = null;
		_parent_pop_sizes = null;
		_age = NO_AGE;
	}

	/**
	 * This constructor instantiates an isolated node. It can also intialize
	 * the node name and data stored in it.
	 *
	 * @param name: The name of the node to be created.
	 * @param data: Data stored in this node.
	 */
	public BniNetNode(String name, T data)
	{
		this();
		setName(name);
		setData(data);
	}

	/**
	 * This function returns the name of the node.
	 */
	public String getName()
	{
		return _name;
	}

	/**
	 * This function returns the age of the node.
	 */
	public double getAge()
	{
		return _age;
	}

	/**
	 * This function returns data stored in this node.
	 */
	public T getData()
	{
		return _data;
	}

	/**
	 * This function computes the indegree of a node. The indegree is equal to the number of
	 * parents it has.
	 */
	public int getIndeg()
	{
		if (_parents == null) {
			return 0;
		}
		else {
			return _parents.size();
		}
	}

	/**
	 * This function computes the outdegree of a node. The outdegree of a node is equal to the number
	 * of children it has.
	 */
	public int getOutdeg()
	{
		if (_children == null) {
			return 0;
		}
		else {
			return _children.size();
		}
	}

	/**
	 * This function checks if a node is the root of the network it belongs.
	 *
	 * @return <code>true</code> if this node is the root; <code>false</code> if it's not.
	 */
	public boolean isRoot()
	{
		return (getIndeg() == 0);
	}


	/**
	 * This function checks if a node is a leaf of a network.
	 *
	 * @return <code>true</code> if this node is a leaf; <code>false</code> otherwise.
	 */
	public boolean isLeaf()
	{
		return (getOutdeg() == 0);
	}

	/**
	 * This function checks if a node is a tree node, i.e. it has only one parent.
	 *
	 * @return <code>true</code> if it is a tree node; <code>false</code> otherwise.
	 */
	public boolean isTreeNode()
	{
		return (getIndeg() < 2);
	}

	/**
	 * This function checks if a node is a network node, i.e. it has more than one parent.
	 *
	 * @return <code>true</code> if this node is a network node; <code>false</code> otherwise.
	 */
	public boolean isNetworkNode()
	{
		return (getIndeg() > 1);
	}

	/**
	 * This function gets all (immediate) child nodes of this node.
	 *
	 * @return an iterable list of children of this node. If this node is a leaf, it returns an empty list.
	 */
	public Iterable<NetNode<T>> getChildren()
	{
		if (_children == null) {
			return new LinkedList<NetNode<T>>();	// Empty list of children.
		}
		else {
			return _children;
		}
	}

	/**
	 * This function gets all (immediate) parents of this node.
	 *
	 * @return an iterable list of parents of this node. If this node has no parents, it returns an empty list.
	 */
	public Iterable<NetNode<T>> getParents()
	{
		if (_parents == null) {
			return new LinkedList<NetNode<T>>();	// Empty list of parents.
		}
		else {
			return _parents;
		}
	}


	/**
	 * This function checks whether the node has a certain parent
	 * */
    public boolean hasParent(NetNode<T> parent){
        if(this.isRoot()){
            return false;
        }
        return _parents.contains(parent);
    }


	/**
	 * This function returns the number of parent of a node.
	 */
	public int getParentCount(){
	    // kliu - this should be guarded appropriately
	    if (_parents == null) {
		    return (0);
	    }
	    else {
		    return (_parents.size());
	    }
	}


    /**
     * This function returns the number of children of a node.
     */
    public int getChildCount(){
        // kliu - this should be guarded appropriately
        if (_children == null) {
            return (0);
        }
        else {
            return (_children.size());
        }
    }


	/**
	 * This function returns the distance from this node to one of its parent.
	 *
	 * @param parent: The parent node that we want to compute the distance from this node to it.

	 * @return the distance from this node to <code>parent</code>. If <code>parent</code> is actually
	 * not a parent of this node, this function returns the constanct <code>Double.NaN</code>.
	 */
	public double getParentDistance(NetNode<T> parent)
	{
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is actually not a parent of this node.
			return Double.NaN;
		}
		//added this for phylonet-hmm, if using pop size and gen time and ages only, then derive this number
		else if(_parent_distances.get(i) != NO_DISTANCE && _parent_distances.get(i).doubleValue() < 0.0){
			System.out.println("Warning: This is not robust! If not using phylonet-hmm and get negative distance by some mistakes, you cannot detect it!");
			return (((BniNetNode)parent).getAge() - _age) / (2.0 * _parent_gen_times.get(i).doubleValue() * _parent_pop_sizes.get(i).doubleValue());
		}
		else {
			return _parent_distances.get(i).doubleValue();
		}
	}

	/**
	 * This function returns the population size from this node to one of its parent.
	 *
	 * @param parent: The parent node that we want to compute the distance from this node to it.

	 * @return the distance from this node to <code>parent</code>. If <code>parent</code> is actually
	 * not a parent of this node, this function returns the constanct <code>Double.NaN</code>.
	 */
	public double getParentPopSize(NetNode<T> parent)
	{
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is actually not a parent of this node.
			return Double.NaN;
		}
		else {
			return _parent_pop_sizes.get(i).doubleValue();
		}
	}

	/**
	 * This function returns the generation time from this node to one of its parent.
	 *
	 * @param parent: The parent node that we want to compute the distance from this node to it.

	 * @return the distance from this node to <code>parent</code>. If <code>parent</code> is actually
	 * not a parent of this node, this function returns the constanct <code>Double.NaN</code>.
	 */
	public double getParentGenTime(NetNode<T> parent)
	{
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is actually not a parent of this node.
			return Double.NaN;
		}
		else {
			return _parent_gen_times.get(i).doubleValue();
		}
	}

	/**
	 * This function updates the data stored in this node.
	 *
	 * @param data: The reference to the new data that we want to store in this node.
	 */
	public void setData(T data)
	{
		_data = data;
	}



	/**
	 * This function changes the name of this node.
	 *
	 * @param name: The new name for this node.
	 */
	public void setName(String name)
	{
		_name = name;
	}



	/**
	 * This function changes the age of this node.
	 *
	 * @param age: The new name for this node.
	 */
	public void setAge(double age)
	{
		_age = age;
	}

	/**
	 * This functions sets the distance from this calling node to <code>parent</code> with
	 * the new value <code>newDistance</code>.
	 *
	 * @param parent: A parent of this node that it wants to modify the distance.
	 * @param distance: New value for the distance from this node to <code>parent</code>.
	 *
	 * @return: <code>true</code> if the operation parent is indeed a parent of this node;
	 * <code>false</code> otherwise.
	 */
	public boolean setParentDistance(NetNode<T> parent, double distance)
	{
		int i = _parents.indexOf(parent);

		if(distance < 0 && distance != NetNode.NO_DISTANCE) {
			//System.out.println(parent + " " + distance + "!!!!!!");
		}

		if (i == -1) {	// parent is not a parent of this node. Operation failed.
			return false;
		}
		else {
			// Update the distance.
			_parent_distances.remove(i);
			_parent_distances.add(i, new Double(distance));

			return true;
		}
	}



	/**
	 * This functions sets the population size from this calling node to <code>parent</code> with
	 * the new value <code>newDistance</code>.
	 *
	 * @param parent: A parent of this node that it wants to modify the distance.
	 * @param popSize: New value for the distance from this node to <code>parent</code>.
	 *
	 * @return: <code>true</code> if the operation parent is indeed a parent of this node;
	 * <code>false</code> otherwise.
	 */
	public boolean setParentPopSize(NetNode<T> parent, double popSize)
	{
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is not a parent of this node. Operation failed.
			return false;
		}
		else {
			// Update the distance.
			if(_parent_pop_sizes != null) {
				_parent_pop_sizes.remove(i);
			}else {
				_parent_pop_sizes = new LinkedList<Double>();
			}
			_parent_pop_sizes.add(i, new Double(popSize));

			return true;
		}
	}

	/**
	 * This functions sets the generation time from this calling node to <code>parent</code> with
	 * the new value <code>newDistance</code>.
	 *
	 * @param parent: A parent of this node that it wants to modify the distance.
	 * @param genTime: New value for the distance from this node to <code>parent</code>.
	 *
	 * @return: <code>true</code> if the operation parent is indeed a parent of this node;
	 * <code>false</code> otherwise.
	 */
	public boolean setParentGenTime(NetNode<T> parent, double genTime)
	{
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is not a parent of this node. Operation failed.
			return false;
		}
		else {
			// Update the distance.
			if(_parent_gen_times != null) {
				_parent_gen_times.remove(i);
			}else
			{
				_parent_gen_times = new LinkedList<Double>();
			}
			_parent_gen_times.add(i, new Double(genTime));

			return true;
		}
	}

	/**
	 * This function connects an existing node (the node that makes a call to this
	 * function) to another node <code>child</code>. The calling code will add <code>child</code>
	 * to its list of children if <code>child</code> has not been already a child of the calling node.
	 *
	 * @param child: The node that the calling node wants to connects to.
	 * @param distance: The distance between the calling code and <code>child</code>.
	 *
	 * @return <code>true</code> if this function succeeded; <code>false</code> otherwise.
	 */
	public boolean adoptChild(NetNode<T> child, double distance)
	{
		assert(child instanceof BniNetNode);


		if (_children == null) {
			_children = new LinkedList<NetNode<T>>();
		}


		BniNetNode<T> ref = (BniNetNode<T>) child;

		if (!_children.contains(ref)) {	// child is not already in this node's list of children.
			// Update info. for this node.
			_children.add(ref);

			// Update info. for child.
			if (ref._parents == null) {
				ref._parents = new LinkedList<NetNode<T>>();
				ref._parent_distances = new LinkedList<Double>();
                ref._parent_support = new LinkedList<Double>();
                ref._parent_probabiliites = new LinkedList<Double>();
			}
			ref._parents.add(this);
			ref._parent_distances.add(new Double(distance));
            ref._parent_support.add(NetNode.NO_SUPPORT);
            ref._parent_probabiliites.add(NetNode.NO_PROBABILITY);

			return true;
		}
		else {	// child is already a child of this node; exit.
			return false;
		}
	}

	/**
	 * This function makes <code>child</code> no longer a child of this node.
	 *
	 * @param child: The node to be removed.
	 *
	 * @return: <code>true</code> if <code>child</code> is indeed a child of this node.
	 */
	public boolean removeChild(NetNode<T> child)
	{
		assert(child instanceof BniNetNode);

		// Delete link from this node to child.
		int index = -1;
		if (_children == null){
			return false;
		}
		index = _children.indexOf(child);
		if(index == -1){
			return false;	// child is not in the list _children; exit.
		}
		_children.remove(index);

		BniNetNode<T> ref = (BniNetNode<T>) child;
		int i = ref._parents.indexOf(this);

		ref._parents.remove(i);
		ref._parent_distances.remove(i);
        ref._parent_probabiliites.remove(i);
        ref._parent_support.remove(i);
		return true;
	}


	/**
	 * This function makes <code>parent</code> no longer a parent of this node.
	 *
	 * @param parent: The node to be removed.
	 *
	 * @return: <code>true</code> if <code>parent</code> is indeed a parent of this node.
	 */
    public boolean removeParent(NetNode<T> parent){
        assert(parent instanceof BniNetNode);

        // Delete link from this node to child.
        if (_parents == null){
            return false;
        }
        int index = _parents.indexOf(parent);
        if(index == -1){
            return false;	// child is not in the list _children; exit.
        }
        _parents.remove(index);
        _parent_distances.remove(index);
        _parent_probabiliites.remove(index);
        _parent_support.remove(index);

        if(_parents.size() == 0){
            _parents = null;
            _parent_distances = null;
            _parent_probabiliites = null;
            _parent_support = null;
        }
        else if(_parents.size()==1 && _children.size()==1){
            NetNode<T> anotherParent = _parents.get(0);
            NetNode<T> child = _children.get(0);
            anotherParent.adoptChild(child, child.getParentDistance(this)+this.getParentDistance(anotherParent));
            child.setParentProbability(anotherParent, this.getParentProbability(anotherParent) * child.getParentProbability(this));
            child.setParentSupport(anotherParent, this.getParentSupport(anotherParent) * child.getParentSupport(this));
            anotherParent.removeChild(this);
            this.removeChild(child);
        }

        BniNetNode<T> ref = (BniNetNode<T>) parent;
        ref._children.remove(this);
        return true;
    }


	/**
	 * This function removes the node itself
	 */
    public void removeItself()
    {
        if(_children!=null){
            List<NetNode<T>> childrenList = new ArrayList<NetNode<T>>();
            for(NetNode<T> child: this.getChildren()){
                childrenList.add(child);
            }
            for(NetNode<T> child: childrenList){
                ((BniNetNode)child).removeParent(this);
            }
        }
        if(_parents!=null){
            List<NetNode<T>> parentList = new ArrayList<NetNode<T>>();
            for(NetNode<T> parent: this.getParents()){
                parentList.add(parent);
            }
            for(NetNode<T> parent: parentList){
                ((BniNetNode)parent).removeChild(this);
            }
        }
    }


	/**
	 * This functions sets the inheritance probability from this calling node to <code>parent</code> with
	 * the new value <code>probability</code>.
	 *
	 * @param parent: 		A parent of this node that it wants to modify the distance.
	 * @param probability: 	New value for the inheritance probability from this node to <code>parent</code>.
	 */
	public void setParentProbability(NetNode<T> parent, double probability)
	{

        if(probability != NetNode.NO_PROBABILITY && (probability < 0 || probability > 1))
        {
            throw new IllegalArgumentException("Probability values must be between zer and one.  Found: " + probability);
        }
		int i = _parents.indexOf(parent);
        
		if (i == -1) {	// parent is not a parent of this node. Operation failed.

            throw new IllegalArgumentException("Passed parent is not a parent of the node.");
		}
		else {

			_parent_probabiliites.set(i, probability);
		}
	}


	/**
	 * This function returns the inheritance probability from this node to one of its parent.
	 *
	 * @param parent: The parent node that we want to compute the inheritance probability from this node to it.

	 * @return the inheritance probability from this node to <code>parent</code>.
	 */
    public double getParentProbability(NetNode<T> parent)
    {
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is not a parent of this node. Operation failed.
            throw new IllegalArgumentException("Passed parent is not a parent of the node.");
		}
		else {
			return _parent_probabiliites.get(i);
		}
	}


	/**
	 * This functions sets the support from this calling node to <code>parent</code> with
	 * the new value <code>probability</code>.
	 *
	 * @param parent: 		A parent of this node that it wants to modify the distance.
	 * @param support: 	New value for the support from this node to <code>parent</code>.
	 */
    public void setParentSupport(NetNode<T> parent, double support)
	{
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is not a parent of this node. Operation failed.

            throw new IllegalArgumentException("Passed parent is not a parent of the node.");
		}
		else {
			// Update the distance.
			_parent_support.remove(i);
			_parent_support.add(i, new Double(support));
		}
	}


	/**
	 * This function returns the support from this node to one of its parent <code>parent</code>.
	 *
	 * @param parent: The parent node that we want to compute the support from this node to it.

	 * @return the support from this node to <code>parent</code>.
	 */
    public double getParentSupport(NetNode<T> parent)
    {
		int i = _parents.indexOf(parent);

		if (i == -1) {	// parent is not a parent of this node. Operation failed.

            throw new IllegalArgumentException("Passed parent is not a parent of the node.");
		}
		else {
			return _parent_support.get(i);
		}
	}

	/**
	 * Get root pop size
	 * @return root pop size
     */
	public double getRootPopSize() {
		return isRoot() ? _root_popSize : NO_ROOT_POPSIZE;
	}

	/**
	 * Set root pop size
	 * @param popSize   root pop size
     */
	public void setRootPopSize(double popSize) {
		if (isRoot()) {
			_root_popSize = popSize;
		} else {
			throw new RuntimeException("Don't allow setting pop size for a node that is not root");
		}
	}


	// Data members
	private T _data;		// Additional data we want to store in this node
	private String _name;	// Node's name
	private List<NetNode<T>> _children;		// List of this node's children
	private List<NetNode<T>> _parents;		// List of parents of this node.
	private List<Double> _parent_probabiliites;           // List of gammas to its children
    public List<Double> _parent_support;           // List of gammas to its children // TODO
	private List<Double> _parent_distances;	// List of distances from this node to its parents.
	//LEO: adding new parameters aug. 9, 2016
	private List<Double> _parent_pop_sizes; //population sizes from this node to its parents
	private List<Double> _parent_gen_times; //generation times from this node to its parents
	private double _age; //the time in years of the node
	// dw20: adding new parameters 08/17/2016
	private double _root_popSize = NO_ROOT_POPSIZE;

}

