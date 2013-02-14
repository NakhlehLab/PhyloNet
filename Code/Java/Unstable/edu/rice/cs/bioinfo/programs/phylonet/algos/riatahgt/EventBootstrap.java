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

package edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt;

import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * This class is to compute bootstrap of HGT events. The class assumes that the bootstrap values are stored
 * in the data field of nodes.
 *
 * @author Cuong Than
 */
public class EventBootstrap {
	// Data members
	private STITree<Double> _st, _gt;							// Species and gene trees.
	private STITree<List<HgtScenario>> _solution;	// Solution tree computed by RIATA.
	private List<HgtEvent> _events;		// Events for reconciling _st and _gt.
	private List<Double> _bootstraps;	// Bootstraps value for events.
	private STITree<Double> _st_prime, _gt_prime;	// Used for computing bootstraps.

	// Member functions
	public EventBootstrap(STITree<Double> st, STITree<Double> gt, STITree<List<HgtScenario>> solution) {
		_st = st;
		_gt = gt;
		_solution = solution;
		_events = new LinkedList<HgtEvent>();
		_bootstraps = new LinkedList<Double>();
	}

	/**
	 * Compute bootstrap values for events in the solution tree.
	 */
	public void computeBootstrap() {
		List<HgtScenario> scenarios = generateScenarios();



		try {
			for (HgtScenario hs : scenarios) {
				computeBootstrapHelper(hs);
			}
		}
		catch (Exception e) {
			System.err.println(e.getMessage());
			return;
		}
	}

	/**
	 * Return all *different* events found by Riata for (_st, _gt).
	 */
	public List<HgtEvent> getEvents() {
		return _events;
	}

	/**
	 * Return bootstrap values for HGT events.
	 */
	public List<Double> getBootstraps() {
		return _bootstraps;
	}

	/**
	 * Generate all possible HGT scenarios.
	 */
	private List<HgtScenario> generateScenarios() {
		// Get all non-empty components.
		List<STINode<List<HgtScenario>>> components = new LinkedList<STINode<List<HgtScenario>>>();
		for (STINode<List<HgtScenario>> node : _solution.getNodes()) {
			if (!node.getData().isEmpty()) {
				components.add(node);
			}
		}

		// Generate all different combinations of HgtScenario from all components.
		List<HgtScenario> res = new LinkedList<HgtScenario>();
		List<Integer> counters = new LinkedList<Integer>();

		if (components.isEmpty()) {
			return res;
		}

		else {
			for (int i = 0; i < components.size(); i++) {
				counters.add(0);
			}

			boolean stop = false;
			while (!stop) {
				// Get one scenario.
				HgtScenario scenario = new HgtScenario();
				for (int i = 0; i < components.size(); i++) {
					int countVal = counters.get(i);
					scenario.addEvents(components.get(i).getData().get(countVal));
				}

				res.add(scenario);

				// Go to the next one.
				int i = components.size() - 1;
				while (i >= 0) {
					int countVal = counters.get(i);
					if (countVal < components.get(i).getData().size() - 1) {
						counters.remove(i);
						counters.add(i, countVal + 1);
						break;
					}
					else {
						counters.remove(i);
						counters.add(i, 0);
						i--;
					}
				}

				if (i < 0) {
					stop = true;
				}
			}

			return res;
		}
	}

	/**
	 * Compute bootstraps for events in this scenario. Update _events and _bootstraps
	 * if there's any new event in <code>scenario</code>. The implementation here is based
	 * on the algorithm in the paper "SPR-...".
	 *
	 * @param scenario: A set of HGT events for reconciling _st and _gt.
	 */

	/**
	 * Compute bootstraps for events in this scenario. Update _events and _bootstraps
	 * if there's any new event in <code>scenario</code>. The implementation here is based
	 * on the algorithm in the paper "SPR-...".
	 *
	 * @param scenario: A set of HGT events for reconciling _st and _gt.
	 */
	private void computeBootstrapHelper(HgtScenario scenario) throws Exception {
		List<TNode> newLeaves = createTempTrees(scenario);



		for (HgtEvent event : scenario.getEvents()) {
			if (event.isBad()) {
				continue;
			}



            // Compute bootstrap by using GT. First, compute the set P.

            STINode<Double> gtDest = _gt_prime.getNode(event.getDestEdge().getName());;
            STINode<Double> gtDestParent =  gtDestParent = gtDest.getParent();
            Set<STINode<Double>> setP = getSubleaves(gtDestParent);





			for (TNode node : newLeaves) {
				if (setP.contains(node)) {
					setP.remove(node);	// Remove unwanted leaves resulted in _gt_prime, which are actually original internal nodes.
				}
			}
			if (setP.isEmpty()) {
				System.err.println("Error when computing the bootstral value: Set P is empty.");
				throw new Exception("Error when computing the bootstral value: Set P is empty.");
			}

			// Compute the set Q.
			Set<STINode<Double>> setQ;
			STINode<Double> stDest = _st_prime.getNode(event.getDestEdge().getName());

			while (true) {	// Repeat until Q is not empty.
				STINode<Double> yPrime = stDest.getParent();	// Y' in ST.
				STINode<Double> gtParent = _gt_prime.getNode(yPrime.getParent().getName());

				setQ = getSubleaves(gtParent);
				for (TNode node : newLeaves) {
					if (setQ.contains(node)) {
						setQ.remove(node);
					}
				}
				if (!setQ.isEmpty()) {
					break;
				}
				else {
					stDest = yPrime;
				}
			}

			if (setQ.isEmpty()) {
				System.err.println("Error when computing the bootstrap value: Set Q is empty.");
				throw new Exception("Error when computing the bootstrap value: Set Q is empty.");
			}

			// Get the path between LCA(P) and LCA(Q), and compute bootstrap value.
			SchieberVishkinLCA lca = new SchieberVishkinLCA(_gt);
			Set<STINode<Double>> temp = new HashSet<STINode<Double>>();

			// p = LCA(P).
			for (TNode node : setP) {
				if (_gt.getNode(node.getName()) != null) {
					temp.add(_gt.getNode(node.getName()));
				}
				else {
					System.err.println("Error when computing the bootstrap value: Leaf " + node.getName() + " does not exist in the gene tree.");
					throw new Exception("Error when computing the bootstrap value: Leaf " + node.getName() + " does not exist in the gene tree.");
				}
			}

			STINode<Double> p = (STINode<Double>) lca.getLCA(temp);
			if (p == null) {
				System.err.println("Error when computing the bootstrap value: Cannot get the common ancestor of P");
				throw new Exception("Error when computing the bootstrap value: Cannot get the common ancestor of P");
			}

			// q = LCA(Q).
			temp.clear();
			for (TNode node : setQ) {
				if (_gt.getNode(node.getName()) != null) {
					temp.add(_gt.getNode(node.getName()));
				}
				else {
					System.err.println("Error when computing the bootstrap value: Leaf " + node.getName() + " does not exist in the gene tree.");
					throw new Exception("Error when computing the bootstrap value: Leaf " + node.getName() + " does not exist in the gene tree.");
				}
			}

			STINode<Double> q = (STINode<Double>) lca.getLCA(temp);
			if (q == null) {
				System.err.println("Error when computing the bootstrap value: Cannot get the common ancestor of Q");
				throw new Exception("Error when computing the bootstrap value: Cannot get the common ancestor of Q");
			}

			// Get the bootstrap value.
			double gtBootstrap = getEdgeBootstrap(p, q);

			// Compute bootstrap in ST.
			p = _st.getNode(event.getSourceEdge().getName());
			q = _st.getNode(event.getDestEdge().getName());

			double stBootstrap;

			if (p == null || q == null) {
				System.err.println("Error when computing the bootstrap value: Cannot get the source and destination nodes in the species tree.");
				throw new Exception("Error when computing the bootstrap value: Cannot get the source and destination nodes in the species tree.");
			}
			else {
				stBootstrap = getEdgeBootstrap(p, q);
			}

			// Update _events and _bootstraps.
			double bs;

			if (stBootstrap <= 0) {
				bs = gtBootstrap;
			}
			else {
				bs = (gtBootstrap < stBootstrap) ? gtBootstrap : stBootstrap;
			}

			if (!_events.contains(event)) {
				_events.add(event);
				_bootstraps.add(bs);
			}
			else {
				int i = _events.indexOf(event);
				if (bs < _bootstraps.get(i)) {
					_bootstraps.remove(i);
					_bootstraps.add(i, bs);
				}
			}
		}
	}

	private void computeBootstrapHelper2(HgtScenario scenario) {
		List<TNode> newLeaves = createTempTrees(scenario);




		for (HgtEvent event : scenario.getEvents()) {
			if (event.isBad()) {
				continue;
			}

			// Compute bootstrap by using GT. First, compute the set P.
			STINode<Double> gtDest = _gt_prime.getNode(event.getDestEdge().getName());
			Set<STINode<Double>> setP = getSubleaves(gtDest.getParent());

			// Compute the set Q.
			STINode<Double> source = _gt_prime.getNode(event.getSourceEdge().getName());
			Set<STINode<Double>> setX = getSubleaves(source.getParent());
			Set<STINode<Double>> setQ;
			STINode<Double> stDest = _st_prime.getNode(event.getDestEdge().getName());

			while (true) {	// Repeat until Q is not empty.
				STINode<Double> yPrime = stDest.getParent();	// Y' in ST.
				STINode<Double> gtParent = _gt_prime.getNode(yPrime.getParent().getName());

				setQ = getSubleaves(gtParent);
				setQ.removeAll(setX);
				gtDest = _gt_prime.getNode(yPrime.getName());
				setQ.removeAll(getSubleaves(gtDest));

				if (setQ.isEmpty()) {
					stDest = yPrime;
				}
				else {
					break;
				}
			}

			// Get the path between LCA(P) and LCA(Q), and compute bootstrap value.
			SchieberVishkinLCA lca = new SchieberVishkinLCA(_gt);
			Set<STINode<Double>> temp = new HashSet<STINode<Double>>();

			for (TNode node : setP) {
				if (_gt.getNode(node.getName()) != null) {
					temp.add(_gt.getNode(node.getName()));
				}
			}

			STINode<Double> p = temp.isEmpty() == true ? null : ((STINode<Double>) lca.getLCA(temp));
			temp.clear();
			for (TNode node : setQ) {
				if (_gt.getNode(node.getName()) != null) {
					temp.add(_gt.getNode(node.getName()));
				}
			}

			STINode<Double> q = temp.isEmpty() == true ? null : ((STINode<Double>) lca.getLCA(temp));
			double gtBootstrap;

			if (p == null || q == null) {
				gtBootstrap = 0;
			}
			else {
				gtBootstrap = getEdgeBootstrap(p, q);
			}

			// Compute bootstrap in ST.
			p = _st.getNode(event.getSourceEdge().getName());
			q = _st.getNode(event.getDestEdge().getName());

			double stBootstrap = getEdgeBootstrap(p, q);

			// Update _events and _bootstraps.
			double bs;

			if (stBootstrap == 0) {
				bs = gtBootstrap;
			}
			else {
				bs = (gtBootstrap < stBootstrap) ? gtBootstrap : stBootstrap;
			}

			// Take the highest bootstrap value for the event.
			if (!_events.contains(event)) {
				_events.add(event);
				_bootstraps.add(bs);
			}
			else {
				int i = _events.indexOf(event);
				if (bs > _bootstraps.get(i)) {
					_bootstraps.remove(i);
					_bootstraps.add(i, bs);
				}
			}
		}
	}

	private List<TNode> createTempTrees(HgtScenario scenario) {
		_st_prime = new STITree<Double>(_st);
		_gt_prime = new STITree<Double>(_st);




		List<TNode> newLeaves = new LinkedList<TNode>();

		int i = 0;
		for (HgtEvent event : scenario.getEvents()) {
			if (!event.isBad()) {



				// Add source and destination nodes to _st_prime.
				TNode source = _st_prime.getNode(event.getSourceEdge().getName());
				TNode dest = _st_prime.getNode(event.getDestEdge().getName());
				TMutableNode temp;

				if (!source.isRoot()) {
                    temp = ((TMutableNode) source.getParent()).createChildWithUniqueName();
					//temp = ((TMutableNode) source.getParent()).createChild("IS_" + i);
					temp.adoptChild((TMutableNode) source);
				}
				else {
                    temp = ((TMutableNode) source).createChildWithUniqueName();
					//temp = ((TMutableNode) source).createChild("IS_" + i);
					temp.makeRoot();
				}

				//temp = ((TMutableNode) dest.getParent()).createChild("ID_" + i);
                temp = ((TMutableNode) dest.getParent()).createChildWithUniqueName();
				temp.adoptChild((TMutableNode) dest);

				// Add source and destination nodes to _gt_prime.
				source = _gt_prime.getNode(event.getSourceEdge().getName());
				dest = _gt_prime.getNode(event.getDestEdge().getName());


				if (!source.isRoot()) {
					temp = ((TMutableNode) source.getParent()).createChild("IS_" + i);
					temp.adoptChild((TMutableNode) source);
				}
				else {
					temp = ((TMutableNode) source).createChild("IS_" + i);
					temp.makeRoot();
				}


				TMutableNode temp2 = ((TMutableNode) dest.getParent()).createChild("ID_" + i);
				temp2.adoptChild((TMutableNode) dest);

				TNode parent = temp2.getParent();	// This node can possibly become a leaf.
				if (parent.getChildCount() == 1) {
					newLeaves.add(parent);
				}



				temp.adoptChild(temp2);

				// Next event.
				i++;


			}
		}



		return newLeaves;
	}

	/**
	 * Get the set of leaves under <code>node</code>.
	 *
	 * @param node: node to find the set of leaves under it.
	 */
	private Set<STINode<Double>> getSubleaves(STINode<Double> node) {

        LinkedList<STINode<Double>> toExplore = new LinkedList<STINode<Double>>();
        toExplore.add(node);

        HashSet<STINode<Double>> subLeafsAccum = new HashSet<STINode<Double>>();

        while(!toExplore.isEmpty())
        {
            STINode<Double> someNode = toExplore.removeFirst();
            if(someNode.isLeaf())
            {
                subLeafsAccum.add(someNode);
            }
            else
            {
                for (STINode<Double> child : someNode.getChildren()) {

                    toExplore.add(child);
                }
            }
        }

        return subLeafsAccum;
	}

	/**
	 * Get the maximium bootstrap value of edges between two nodes p and q.
	 *
	 * @param p, q
	 */
	private double getEdgeBootstrap(STINode<Double> p, STINode<Double> q) {
		if (p == q) {
			return TNode.NO_DISTANCE;
		}

		Set<STINode<Double>> temp = getSubleaves(p);
		temp.addAll(getSubleaves(q));

		SchieberVishkinLCA lca = new SchieberVishkinLCA(p.getTree());
		STINode<Double> cross = (STINode<Double>) lca.getLCA(temp);

		double max = TNode.NO_DISTANCE;
		while (p != cross) {
			// double d = p.getParentDistance();
			if (p.getData() != null) {
				double d = p.getData();

				if (d != TNode.NO_DISTANCE && d > max) {
					max = d;
				}

				p = p.getParent();
			}
			else {
				break;
			}
		}

		while (q != cross) {
			// double d = q.getParentDistance();
			if (q.getData() != null) {
				double d = q.getData();
				if (d != TNode.NO_DISTANCE && d > max) {
					max = d;
				}

				q = q.getParent();
			}
			else {
				break;
			}
		}

		return max;
	}
}

