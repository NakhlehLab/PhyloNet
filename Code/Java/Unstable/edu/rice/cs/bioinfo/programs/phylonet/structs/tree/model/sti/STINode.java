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

package edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.PostTraversal;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;


public class STINode<D extends Object> implements TMutableNode {
	private D _data;
	protected STITree<D> _tree;
	protected boolean _valid = true;
	protected int _id;
	protected String _name = STITree.NO_NAME;
	protected STINode<D> _parent;
	protected LinkedList<STINode<D>> _children = new LinkedList<STINode<D>>();
	protected double _distance = TMutableNode.NO_DISTANCE;

	
	
	/**
	 * Constructor
	 */
	protected STINode(STITree<D> tree, int id, String name, STINode<D> parent, D data) {
		if (tree != null) {
			while (tree._nodes.containsKey(id)) {
				id++;
			}
		}
	
		_id = id;
		_tree = tree;
		_parent = parent;
		_data = data;
			
		setName(name);	
		
		_tree._nodes.put(getID(), this);
		_tree._node_set.add(this);		
				
		if (_tree._node_set.size() != _tree._nodes.size()) {
			System.err.println("STINode: Inconsistent _node_set and _nodes");
		}
	}


    /**
     * This function computes the depth of this node
     */
	public int getHeight() {
		if(this.isLeaf()) {
			return 1;
		} else {
			int max_height = 0;
			
			for(STINode child : _children) {
				int h = child.getHeight();
				
				if(h > max_height) {
					max_height = h;
				}
			}
			return max_height + 1;
		}
	}


    /**
     * This function returns the siblings of this node
     */
	public List<TNode> getSiblings(){
		List<TNode> siblings = new ArrayList<TNode>();
		for(TNode n :this.getParent().getChildren()){
			if(n.getID()!=this.getID()){
				siblings.add(n);
			}
		}
		return siblings;
	}
    
    
    /**
     * This function returns the data of this node
     */
	public D getData() { return _data; }

    
    /**
     * This function sets the data of this node
     */
	public void setData(D data) { _data = data; }


    /**
     * This function returns if this node is valid or not
     */
	public boolean isValid() {
		return _valid;
	}


    /**
     * This function returns the id of this node
     */
	public int getID() {
		return _id;
	}


    /**
     * This function returns the name of this node
     */
	public String getName() {
		return _name;
	}

    
    /**
     * This function sets the name of this node
     * Duplicate names are not allowed
     */
	public void setName(String name) {
		if(!name.equals(STITree.NO_NAME) && _tree._name2node.containsKey(name)) {
			throw new RuntimeException("Illegal assignment of duplicate node name '" + name + "'");
		}
		
		_tree._name2node.remove(_name);
		
		_name = name;
		
		_tree._name2node.put(_name, this);
	}


    /**
     * This function returns the tree this node belongs to
     */
	public STITree<D> getTree() {
		if(_valid) {
			return _tree;
		} else {
			return null;
		}
	}


    /**
     * This function returns the parent of this node
     */
	public STINode<D> getParent() {
		return _parent;
	}


    /**
     * This function creates a child with no name of this node
    */
	public STINode<D> createChild() {
		return createChild(STITree.NO_NAME);
	}


    /**
     * This function creates a child with unique name of this node
     */
    public STINode<D> createChildWithUniqueName()
    {
        String prefix = "n";

        for(int i = 0; true; i++)
        {
            String name = prefix + i;

            if(! this._tree._name2node.keySet().contains(name))
            {
                return createChild(name);
            }
        }
    }


    /**
     * This function creates a child with name <code>name<code/>
     */
	public STINode<D> createChild(String name) {
		
		STINode<D> child = new STINode<D>(_tree, _tree._next_node_id++, name, this, null);
		
		_children.add(child);
		
		_tree._leaf_count = STITree.UNCOUNTED;
		
		return child;
	}


    /**
     * This function adopts a child with <code>nchild<code/>
     */
	public void adoptChild(TMutableNode nchild) {
		STINode<D> child = (STINode<D>) nchild;

        if(child._tree == null)
        {
            for(TNode c : nchild.postTraverse())
            {
               STINode mc = (STINode)c;
               mc._tree = _tree;
               _tree._nodes.put(mc._id, mc);
               _tree._node_set.add(mc);
               if(mc._name != null)
                   _tree._name2node.put(mc._name, mc);
            }
        }

		// sanity check to make sure that the child is a member of this tree
		assert child._tree == _tree;
        
		if(child._parent == this) {
			return;
		}
		//child is not root
		if(child._parent != null) {			
			child._parent._children.remove(child);
		}
		
		child._parent = this;
		_children.add(child);
		
		_tree._leaf_count = STITree.UNCOUNTED;
	}


    /**
     * This function sets the parent <code>newParent<code/> of this node
     */
    public void setParent(TMutableNode newParent)
    {
        if(isRoot())
        {
            throw new IllegalStateException("Roots may not change parents.");
        }

        this.getParent().removeChild(this, false);
        _parent = null;
        newParent.adoptChild(this);
    }


    /**
     * This function removes edge <code>adjacentNode<code/> adjacent to this node
     */
    public void removeEdge(TMutableNode adjacentNode)
    {
        if(adjacentNode.equals(_parent))
        {
            _parent.removeChild(this, false);
            this._tree = null;

            for(TNode d : this.postTraverse())
            {
                ((STINode)d)._tree = null;
            }
        }
        else
        {
            removeChild(adjacentNode, false);
            STINode adjST = (STINode)adjacentNode;

            adjST._tree = null;

            for(Object d : adjST.postTraverse())
            {
                ((STINode)d)._tree = null;
            }
        }
    }


    /**
     * This function removes child <code>child<code/> while adopting all its children if <code>adopt_all</code> is set to true
     */
	public void removeChild(TMutableNode child, boolean adopt_all) {
		STINode<D> batichild = (STINode<D>) child;
		
		if(adopt_all == true) {
			// grab all the nodes children and make this node their parent
			Iterator<STINode<D>> i = batichild._children.iterator();
			
			while(i.hasNext()) {
				i.next()._parent = this;
			}
			
			// make them our children
            for(STINode<D> childToAdopt: batichild._children) {
                _children.add(childToAdopt);
                childToAdopt.setParentDistance(childToAdopt.getParentDistance() + batichild.getParentDistance());
            }
			
			// remove them from the former parent
			batichild._children.clear();
		}

		_children.remove(batichild);
		
		batichild.removeSelf();	
		
		_tree._leaf_count = STITree.UNCOUNTED;
	}
    

	/**
	 * This method removes the node and any children it has from the STITree object.
	 */
	protected void removeSelf() {
		
		// remove us as a child
		if(_parent != null) {
			_parent._children.remove(this);
		}
		
		// remove our children
		if(!isLeaf()) {
			Iterator i = new LinkedList<STINode<D>>(_children).iterator();
			
			while(i.hasNext()) {
				((STINode) i.next()).removeSelf();
			}
		}

		_tree.removeNodeRecord(this);
		_tree._leaf_count = STITree.UNCOUNTED;
		_valid = false;
	}


    /**
     * This method removes the node and any children it has from the STITree object.
     */
    private void removeSelf2(boolean iterator) {
        // remove us as a child
        if(_parent != null && iterator) {
            _parent._children.remove(this);
        }

        // remove our children
        if(!isLeaf()) {
            Iterator i = new LinkedList<STINode<D>>(_children).iterator();

            while(i.hasNext()) {
                ((STINode) i.next()).removeSelf2(false);
            }
        }

        _tree.removeNodeRecord(this);

        _tree._leaf_count = STITree.UNCOUNTED;
        _valid = false;
    }

    /**
     * This method returns all children of this node
     */
	public Iterable<STINode<D>> getChildren() {
		return new Iterable<STINode<D>>() {
			public Iterator<STINode<D>> iterator() {
				return _children.iterator();
			}
		};
	}


    /**
     * This method returns all nodes adjacent to this node
     */
    public Iterable<STINode<D>> getAdjacentNodes()
    {
        if(isRoot())
        {
            return getChildren();
        }
        else
        {
            LinkedList<STINode<D>> accum = new LinkedList<STINode<D>>();
            for(STINode<D> node : getChildren())
            {
                accum.add(node);
            }

            accum.add(getParent());

            return accum;
        }
    }


    /**
     * This method returns the number of nodes adjacent to this node
     */
    public int getAdjacentNodeCount()
    {
        if(isRoot())
        {
            return getChildCount();
        }
        else
        {
            return getChildCount() + 1;
        }
    }


    /**
     * This method removes all children of this node
     */
	public List<STINode<D>> removeAllChildren(){
		List<STINode<D>> children = new ArrayList<STINode<D>>(_children);
		for(STINode<D> child: children){			
			_children.remove(child);		
			child.removeSelf2(true);
			_tree._leaf_count = STITree.UNCOUNTED;
		}
		return children;
	}


    /**
     * This method get the number of children of this node
     */
	public int getChildCount() {
		return _children.size();
	}


    /**
     * This method gets the number of degrees of this node (indegree + outdegree)
     */
    public int getDegree()
    {
        if(isRoot())
        {
            return _children.size();
        }
        else
        {
            return _children.size() + 1;
        }
    }


    /**
     * This method returns whether this node is leaf or not
     */
	public boolean isLeaf() {
		return _children.isEmpty();
	}


    /**
     * This method returns whether this node is root or not
     */
	public boolean isRoot() {
		return (_parent == null);
	}


    /**
     * This method makes this node the root of the tree
     */
	public void makeRoot() {
		if(this == _tree._root) {
			return;
		}
		
		STINode<D> prev = null;
		STINode<D> node = this;
		STINode<D> next;
		
		// walk up to the root and reverse relationships until we reach the current root
		while(node != null) {
			next = node._parent;
			
			node._parent = prev;
			
			if(next != null) {
				node._children.add(next);
				next._children.remove(node);
			}
			
			prev = node;
			node = next;
		}
		
		// does the old root need to be removed?
		/*
		 * This removal code doesn't need to be here.  An unrooted tree can be appropriately
		 * modeled by choosing an arbitary node in the tree to be the root.  This code
		 * accomodated for a tree being unrooted and specifying an extra 'root' node.  
		 * The user of this class should handle this.
		if(remove_binary_nodes) {
			prev.removeChild(node, true);
		}
		*/
		
		// set this node as the new root
		_tree._root = this;
		_tree._leaf_count = STITree.UNCOUNTED;
	}


    /**
     * This method returns the length of the branch incident with this node
     */
	public double getParentDistance() {
		return _distance;
	}


    /**
     * This method sets the length of the branch incident with this node
     */
	public void setParentDistance(double distance) {
		_distance = distance;
	}


    /**
     * This method returns the the number of leaves under this node
     */
	public int getLeafCount() {
		
		if(isLeaf()) {
			return 1;
		} else {
		
			int num_leaves = 0;
			
			Iterator i = _children.iterator();
			while(i.hasNext()) {
				num_leaves += ((STINode) i.next()).getLeafCount();
			}
			
			return num_leaves;
		}
	}


    /**
     * This method removes the node itself
     */
	public void removeNode() {
		removeSelf();
	}


    /**
     * This method returns the newick string of this node
     */
	public String toNewick() {
		synchronized(STITree.SWRITER) {
			STITree.SWRITER.reset();
			//STITree.NWRITER.writeTree(this, false);
			STITree.NWRITER.writeTree(this, false);
		
			return STITree.SWRITER.toString();
		}
	}


    /**
     * This method creates a child node <code>clade</code>
     */
	public TMutableNode createChild(TNode clade) {
		STINode<D> node = createChild(clade.getName());
        node.setParentDistance(clade.getParentDistance());
		if(((STINode<D>)clade).getData()!=null){
			node.setData(((STINode<D>)clade).getData());
		}
		
		for(TNode child : clade.getChildren()) {
			node.createChild(child);
		}
		
		return node;
	}


	/**
	 * The string version of a node is its name, if it has one.  If it doesn't have a name, then
	 * it prints the subtree it defines.
	 */
	public String toString() { 
		if(getName() != STITree.NO_NAME) {
			return getName();
		} else {
			return toNewick();
		}
	}


    /**
     * This method returns the newick string of this node
     */
	public String toString(int format) {
		switch(format) {
		case Tree.NEWICK_FORMAT:
			return toNewick();
		default:
			throw new RuntimeException("Unknown format " + format);
		}
	}


    /**
     * This method returns whether this node is an ancestor of n
     */
	public boolean isAncestor(TNode n) {
		while(n != null) {
			if(n == this) {
				return true;
			}
			
			n = n.getParent();
		}
		return false;
	}


    /**
     * This method returns all the leaves under this node
     */
	public Iterable<TNode> getLeaves() {
		PostTraversal<D> traversal = new PostTraversal<D>(this);
		List<TNode> leaves = new LinkedList<TNode>();
		
		for (TNode node : traversal) {
			if (node.isLeaf()) {
				leaves.add(node);
			}
		}
		
		return leaves;
	}


    /**
     * This method returns nodes in the post traversal order
     */
	public Iterable<TNode> postTraverse() {
		return new PostTraversal<D>(this);
	}



    /*
    Sort the child list. Use it with caution! Make sure than the orderedChildren contains the same children
    */
    public boolean sortChildren(List<STINode<D>> orderedChildren){
        if(orderedChildren.size()!=_children.size()){
            return false;
        }
        _children.clear();
        for(STINode<D> child: orderedChildren){
            _children.add(child);
        }
        return true;
    }
}
