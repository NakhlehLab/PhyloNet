package edu.rice.cs.bioinfo.programs.phylonet.algos.gtdistribution;


import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.math.BigInteger;
import java.util.*;


public class GTDistributionInNetwork {
	boolean print = false;
	boolean init = false;
	List<String> net_taxa;
	List<String> st_taxa;
	private boolean [][]  R, M, S;
	Map<String,Integer> nname2tamount;  //map the node name in the network to the number of corresponding nodes in the species tree
	Map<String,String> tname2nname;	 //map the node name in the species tree to the name of the corresponding node in the network
	Map<String,List<TNode>> hname2tnodes;  //map the name of hybrid node to the corresponding nodes in the species tree
	Tree st;
	
	/**
	 * Constructor that initialize the variables.
	 */	
	public GTDistributionInNetwork(){
		net_taxa = new ArrayList<String>();
		st_taxa = new ArrayList<String>();
		nname2tamount = new HashMap<String,Integer>(); 
		hname2tnodes = new HashMap<String,List<TNode>>();
		tname2nname = new HashMap<String,String>();
	}
	
	/**
	 * The public function for calculating the probabilities.
	 * @param	net 	the given network
	 * @param 	gts		the given set of gene trees
	 * @param	alleles2species		the mapping from the names of allels to the names of the species. It is used for multiple alleles
	 * @return	a list of probabilities corresponding to the list of gene trees.
	 */
	public List<Double> calculateGTDistribution(Network<Double> net, List<Tree> gts, Map<String,String> allele2species){
		networkToTree(net);
		
		for(NetNode leaf: net.getLeaves()){
			net_taxa.add(leaf.getName());
		}
		for(Map.Entry<String, Integer> entry: nname2tamount.entrySet()){
			for(int i=1; i<=entry.getValue(); i++){
				String name = entry.getKey();
				tname2nname.put(name+"_"+i, name);
			}
		}
		S = calculateSorR(st);
		computeNodesUnderHybrid(st);
		
		List<Double> probs = performCalculating(gts, allele2species);
		return probs;
	}
	
	
	
	/**
	 * The actual calculation of the gene tree probabilities in the network
	 * @param 	gts		the given set of gene trees
	 * @param	alleles2species		the mapping from the names of allels to the names of the species. It is used for multiple alleles
	 * @return	a list of probabilities corresponding to the list of gene trees.
	 */
	
	private List<Double> performCalculating(List<Tree> gts, Map<String,String> allele2species){		
		List<Double> problist = new ArrayList<Double>();
		double totalprob = 0;
		for(Tree gt: gts){		
			System.out.print(gt.toNewick()+": ");
			if(print){
				System.out.println();
			}
			List<String> gt_taxa = Arrays.asList(gt.getLeaves());
			R = calculateSorR(gt);
			
			List<int[]> allmappings = new ArrayList<int[]>();
			allmappings.add(new int[gt_taxa.size()]);
			for(int i=0; i<gt_taxa.size(); i++){
				String gtleaf = gt_taxa.get(i);
				String hleaf;
				if(allele2species!=null){
					hleaf = allele2species.get(gtleaf);
				}
				else{
					hleaf = gtleaf;
				}
				List<int[]> temp = new ArrayList<int[]>();
				temp.addAll(allmappings);
				allmappings.clear();
				for(int j=1; j<=nname2tamount.get(hleaf); j++){
					int index = st_taxa.indexOf(hleaf+"_"+j);
					for(int[] gtl2stl: temp){
						int[] new_map = gtl2stl.clone();
						new_map[i] = index;
						allmappings.add(new_map);
					}
				}
			}
			
			double gtprob = 0;
			for(int[] mapping: allmappings){
				
				if(print){
					for(int i=0; i<mapping.length; i++){
						System.out.print(gt_taxa.get(i)+"->"+st_taxa.get(mapping[i])+"\t");
					}
					System.out.println();
				}
				
				Map<String,String> aname2tname = new HashMap<String,String>();
				for(int i = 0; i<mapping.length; i++){
					aname2tname.put(gt_taxa.get(i), st_taxa.get(mapping[i]));
				}
				calculateM(gt,st,aname2tname);
				Map<CEPair, Integer> ro = new HashMap<CEPair,Integer>();
				List<int[]> histories = new ArrayList<int[]>();
				int[] his = new int[gt.getNodeCount()];
				Arrays.fill(his, -1);
				histories.add(his);
				enumCoalHistories(gt.getRoot(), ro, histories);
				ro.clear();
				double gtmapprob = 0;
				for(int[] history: histories){
					double gtmaphisprob = 1;
					boolean first = true;
					for(TNode b: st.postTraverse()){						
						String nname = tname2nname.get(b.getName());
						if(hname2tnodes.containsKey(nname)){
							continue;
						}
						int u = calculateU(st, b, mapping, history);
						if(u==0)continue;
						int c = calculateC(b,history);
						double gij = gij(b.getParentDistance(),u,u-c);
						long w = calculateW(b,c,history);
						long d = calculateD(u,c);
						double gamma = ((STINode<Double>)b).getData();
						gtmaphisprob *= gij*w/d*Math.pow(gamma, u-c);
						if(print){
							String prefix = "*";
							if(first){
								prefix = "+";
							}
							if(gij!=1){
								System.out.print(prefix+"g"+u+(u-c)+"("+b.getParentDistance()+")");
								first = false;
							}
							if(d!=1){
								System.out.print(prefix+"("+w+"/"+d+")");
								first = false;
							}
							if(gamma!=1 && u-c!=0){
								if(u-c!=1){
									System.out.print(prefix+"("+gamma+")^"+(u-c));
								}
								else{
									System.out.print(prefix+"("+gamma+")");
								}
								first = false;
							}
						}
					}		
					for(Map.Entry<String,List<TNode>> entry: hname2tnodes.entrySet()){
						int sum_u = 0;
						int sum_c = 0;
						double prod_w =1;
						double distance = 0;
						for(TNode hnode: entry.getValue()){
							distance = hnode.getParentDistance();
							int u = calculateU(st, hnode, mapping, history);
							int c = calculateC(hnode,history);
							double gamma = ((STINode<Double>)hnode).getData();
							double w = calculateHW(hnode,c,history);
							gtmaphisprob *= Math.pow(gamma, u-c);
							if(print){
								String prefix = "*";
								if(first){
									prefix = "+";
								}
								if(gamma!=1 && u-c!=0){
									if(u-c!=1){
										System.out.print(prefix+"("+gamma+")^"+(u-c));
									}
									else{
										System.out.print(prefix+"("+gamma+")");
									}
									first = false;
								}
							}
							sum_u += u;
							sum_c += c;
							prod_w *= w;
						}
						
						double gij = gij(distance,sum_u,sum_u-sum_c);
						long d = calculateD(sum_u,sum_c);						
						prod_w *= fact(1,sum_c);
						gtmaphisprob *= gij*prod_w/d;
						if(print){
							String prefix = "*";
							if(first){
								prefix = "+";
							}
							if(gij!=1){
								System.out.print(prefix+"g"+sum_u+(sum_u-sum_c)+"("+distance+")");
								first = false;
							}
							if(d!=1){
								System.out.print(prefix+"("+prod_w+"/"+d+")");
								first = false;
							}
						}	
					}
					
					if(print)
						System.out.println("");
					
					gtmapprob += gtmaphisprob;
				}
				gtprob += gtmapprob;
				
				if(print)
					System.out.println("");
				
			}
			System.out.println(gtprob);
			problist.add(gtprob);
			totalprob += gtprob;
		}
		System.out.println("total probability:"+totalprob);
		return problist;
		
	}
	
	
	
	/**
	 * The function is to convert a network to a multilabel tree.
	 * @param	net 	the given network
	 */
	private void networkToTree(Network<Double> net){
		adjustNetwork(net);
		st = new STITree<Double>();
		((STINode<Double>)(st.getRoot())).setData(1.0);
		Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();	
		Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
		source.offer(net.getRoot());
		dest.offer((TMutableNode)st.getRoot());
		while(!source.isEmpty()){
			NetNode<Double> parent = source.poll();
			TMutableNode peer = dest.poll();
			
			int index = 0;
			for (NetNode<Double> child : parent.getChildren()) {
				
				TMutableNode copy;
				if (child.getName() == NetNode.NO_NAME) {
					copy = peer.createChild(TNode.NO_NAME);
				}
				else {
					Integer amount = nname2tamount.get(child.getName());
					if(amount==null){
						amount = 0;
					}
					nname2tamount.put(child.getName(), ++amount);
					String newname = child.getName() + "_" + amount;
					copy = peer.createChild(newname);
					if(child.isLeaf()){
						st_taxa.add(newname);
					}
					
					if(child.getParentCount()>1){
						List<TNode> corresNodes = hname2tnodes.get(child.getName());
						if(corresNodes==null){
							corresNodes = new ArrayList<TNode>();
							corresNodes.add(copy);
							hname2tnodes.put(child.getName(), corresNodes);
						}
						else{
							corresNodes.add(copy);
						}
					}
					
				}
				// Update the distance and data for this child.
				double distance = child.getParentDistance(parent);					
				if (distance == NetNode.NO_DISTANCE) {
					//copy.setParentDistance(TNode.NO_DISTANCE);
					copy.setParentDistance(0);
				}
				else {
					copy.setParentDistance(distance);
				}
				
			//	double gamma = parent.getGamma(index);
			//	((STINode<Double>)copy).setData(gamma);
					
				// Continue to iterate over the children of nn and tn.
				source.offer(child);
				dest.offer(copy);
				index ++;
			}
		}
	}
	
	
	/**
	 * The function is to collect all nodes under hybridization so that they can be treated differently when calculating probabilities
	 * @param	st	a tree
	 */
	private void computeNodesUnderHybrid(Tree st){
		Map<Integer, BitSet> map = new HashMap<Integer, BitSet>();
		Map<String, List<TNode>> hname2tnodes2 = new HashMap<String,List<TNode>>();
		for (TNode node : st.postTraverse()) {
			BitSet bs = new BitSet(net_taxa.size());
			if(node.isLeaf()){
				String nname = tname2nname.get(node.getName());
				bs.set(net_taxa.indexOf(nname));
			}
			else{
				for(TNode child : node.getChildren()) {					
					BitSet childBS = map.get(child.getID());
					bs.or(childBS);
				}	
			}
			map.put(node.getID(), bs);
			
			String name = tname2nname.get(node.getName());
			if(name == null){
				name = bs.toString();
			}
			if(hname2tnodes.containsKey(name)){
				continue;
			}
			List<TNode> peerNodes = hname2tnodes2.get(name);
			if(peerNodes!=null){
				peerNodes.add(node);
				if(!node.isLeaf()){
					((STINode<Double>)node).setName(name+"_"+peerNodes.size());
					tname2nname.put(name+"_"+peerNodes.size(), name);
				}
				continue;
			}
			boolean add = false;
			for(Map.Entry<String, List<TNode>> entry: hname2tnodes.entrySet()){
				for(TNode hnode: entry.getValue()){
					if(S[hnode.getID()][node.getID()]){
						add = true;
						break;
					}
				}
			}
			if(add){
				peerNodes = new ArrayList<TNode>();
				peerNodes.add(node);
				if(!node.isLeaf()){
					((STINode<Double>)node).setName(name+"_"+peerNodes.size());
					tname2nname.put(name+"_"+peerNodes.size(), name);
				}
				hname2tnodes2.put(name, peerNodes);
			}
		}
		
		hname2tnodes.putAll(hname2tnodes2);
	}
	
	
	/**
	 * The function is to calculate the g_{ij} function. 
	 * @param	length	the branch length
	 * @param 	i	the number of lineages in
	 * @param	j	the number of lineages out
	 * @return	the resulting probability
	 */	
	private double gij(double length, int i, int j){
		if(length == (TNode.NO_DISTANCE)){
			return 1;
		}
		if(length==0){
			if(i==j)
				return 1;
			else
				return 0;
		}
		if(i==0){
			return 1;
		}
		double result = 0;
		for(int k=j; k<=i; k++){
			result += Math.exp(0.5*k*(1.0-k)*length)*(2.0*k-1.0)*Math.pow(-1,k-j)*fact(j,j+k-2)*fact(i-k+1,i)/(fact(1,j)*fact(1,k-j)*fact(i,i+k-1)); 
		}
		return result;
	}
	
	
	/**
	 * The function is to calculate factorial 
	 * @param	start	the first number
	 * @param 	end		the last number
	 * @return	the resulting factorial
	 */
	private long fact(int start, int end){
		long result = 1;
		for(int i=start; i<=end; i++){
			result = result * i;
		}
		return result;
	}
	

	/**
	 * The function is to calculate "N choose K"
	 */
	private long choose(int N, int K) {
	    BigInteger ret = BigInteger.ONE;
	    for (int k = 0; k < K; k++) {
	        ret = ret.multiply(BigInteger.valueOf(N-k))
	                 .divide(BigInteger.valueOf(k+1));
	    }
	    return ret.longValue();
	}

	
	/**
	 * The function is to calculate the R matrix for the given tree.
	 * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
	 */
	private boolean[][] calculateSorR(Tree tree){
		int nnode = tree.getNodeCount();
		boolean[][] matrix = new boolean[nnode][nnode];
		for(int i=0; i<nnode; i++){
			for(int j=0; j<nnode; j++){
				matrix[i][j] = false;
			}
		}
		Map<Integer, BitSet> map = new HashMap<Integer, BitSet>();
		for (TNode node : tree.postTraverse()) {
			BitSet bs = new BitSet(nnode);
			for(TNode child : node.getChildren()) {
				BitSet childBS = map.get(child.getID());
				
				bs.or(childBS);
			}					
			int id = node.getID();
			for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i+1)) {
			     matrix[id][i] = true;
			 }
			bs.set(id);
			map.put(id, bs);
		}
		return matrix;
	}
	
	
	/**
	 * The function is to calculate the M matrix for the given tree.
	 * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
	 */
	private void calculateM(Tree gt, Tree st, Map<String,String> allele2species){
		int ngtnode = gt.getNodeCount();
		int nstnode = st.getNodeCount();
		M = new boolean[nstnode][ngtnode];
			
		Map<Integer, BitSet> gt_id2bs = new HashMap<Integer, BitSet>();
		List<Integer> rmlist = new ArrayList<Integer>();
		for (TNode node : gt.postTraverse()) {
			BitSet bs = new BitSet(st_taxa.size());
			if (node.isLeaf()) {
				String name = allele2species.get(node.getName());	
				bs.set(st_taxa.indexOf(name));		
				rmlist.add(node.getID());
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = gt_id2bs.get(child.getID());
					bs.or(childCluster);
				}									
			}		
			gt_id2bs.put(node.getID(), bs);
		}
		
		for(Integer rmid: rmlist){
			gt_id2bs.remove(rmid);
		}
		rmlist.clear();
		
		Map<Integer, BitSet> st_id2bs = new HashMap<Integer, BitSet>();
		for (TNode node : st.postTraverse()) {
			BitSet bs = new BitSet(st_taxa.size());
			int st_id = node.getID();
			if (node.isLeaf()) {
				String name = node.getName();			
				bs.set(st_taxa.indexOf(name));						
			}
			else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = st_id2bs.get(child.getID());
					bs.or(childCluster);
				}					
			}
			st_id2bs.put(st_id, bs);	
			
			for(Map.Entry<Integer, BitSet> gt_entry: gt_id2bs.entrySet()){
				int gt_id = gt_entry.getKey();
				BitSet gt_bs = (BitSet) gt_entry.getValue().clone();
				gt_bs.and(bs);
				if(gt_bs.equals(gt_entry.getValue())){
					rmlist.add(gt_id);
					M[st_id][gt_id] = true;
					for(int i=0; i<S.length; i++){
						if(S[i][st_id]){
							M[i][gt_id] = true;
						}
					}
				}
			}
			for(int rmid: rmlist){
				gt_id2bs.remove(rmid);
			}
		}
	}
	
	
	/**
	 * The function is to enumerate coalescent histories.
	 * Read the paper "Gene tree distributions under the coalescent process." by James Degnan to for details.
	 */
	private void enumCoalHistories(TNode gt_cluster, Map<CEPair,Integer> ro, List<int[]> histories){
		if(gt_cluster.isLeaf()) {
			return;
		}
		// if this is the lowest cluster
		if(gt_cluster.getLeafCount() <= 2) {
			List<int[]> temp = new ArrayList<int[]>();
			temp.addAll(histories);
			histories.clear();
			for(int edge=0; edge<M.length; edge++){				
				if(M[edge][gt_cluster.getID()]){
					ro.put(new CEPair(gt_cluster.getID(),edge),1);
					for(int[] his: temp){
						int[] new_his = his.clone();
						new_his[gt_cluster.getID()] = edge;
						histories.add(new_his);
					}
				}
			}
		} else {		
			// compute children first
			for(TNode child : gt_cluster.getChildren()) {
				enumCoalHistories(child, ro, histories);
			}
			
			// compute this cluster's score
			List<int[]> tempHis = new ArrayList<int[]>();
			tempHis.addAll(histories);
			histories.clear();
			
			for(int edge=0; edge<M.length; edge++){
				if(!M[edge][gt_cluster.getID()]){
					continue;
				}
				// compute the product of the sums
				CEPair cp = new CEPair();
				int sum_prod = 1;
				List<int[]> childcoaledges = new ArrayList<int[]>();
				int[] coaledge = new int[R.length];
				Arrays.fill(coaledge, -1);
				childcoaledges.add(coaledge);
				for(TNode child : gt_cluster.getChildren()) {		
					if(sum_prod == 0){
						break;
					}
					// skip leaves since they aren't clusters
					if(child.isLeaf()) {
						continue;
					}
					List<int[]> tempcoaledges = new ArrayList<int[]>();
					tempcoaledges.addAll(childcoaledges);
					childcoaledges.clear();
					
					int sum = 0;
					
					for(int cedge=0; cedge<M.length; cedge++){
						if(M[cedge][child.getID()]){
							if(S[edge][cedge] || edge==cedge){
								cp.set(child.getID(),cedge);
								sum += ro.get(cp);
								for(int[] edges: tempcoaledges){
									int[] newedge = edges.clone();
									newedge[child.getID()] = cedge;
									childcoaledges.add(newedge);
								}
							}
						}			
					}
					tempcoaledges.clear();
					sum_prod *= sum;
				}
				for(int[] history: tempHis){
					for(int[] childcoaledge: childcoaledges){
						boolean add = true;
						for(int i=0; i<childcoaledge.length; i++){
							if(childcoaledge[i]==-1)continue;
							if(history[i]!=childcoaledge[i]){
								add = false;
								break;
							}
						}
						if(add){
							int[] newHis = history.clone();
							newHis[gt_cluster.getID()] = edge;
							histories.add(newHis);
						}
					}
				}
				childcoaledges.clear();
				ro.put(new CEPair(gt_cluster.getID(),edge), Math.max(1,sum_prod));				
			}
			tempHis.clear();
		}
	}
	
	
	/**
	 * The function is to calculate the number of lineages going into a branch
	 * @param	st		the multilabel species tree
	 * @param	node	the node that the branch is incident into
	 * @param	mapping		the mapping
	 * @param	history		the coalescent history of the gene tree
	 */
	private int calculateU(Tree st, TNode node, int[] mapping, int[] history){
		int u = 0;
		for(int i=0; i<mapping.length; i++){
			int mapping_id = st.getNode(st_taxa.get(mapping[i])).getID();
			if(node.isLeaf()){
				if(node.getID() == mapping_id){
					u++;
				}
			}
			else{
				if(S[node.getID()][mapping_id]){
					u++;
				}
			}
		}
		for(int i=0; i<history.length; i++){
			if(history[i]!=-1){
				if(S[node.getID()][history[i]]){
					u--;
				}
			}
		}
		return u;
	}
	
	
	/**
	 * The function is to calculate the number of coalescent events in a branch
	 * @param	node	the node that the branch is incident into
	 * @param	history		the coalescent history
	 */
	private int calculateC(TNode node, int[] history){
		int c = 0;
		for(int i=0; i<history.length; i++){
			if(history[i]==node.getID()){
				c++;
			}
		}
		return c;
	}
	
	
	/**
	 * The function is to calculate the number of possible ordering coalescent events
	 * @param	u	the number of lineages entering the branch
	 * @param	c	the number of coalescent events
	 */
	private long calculateD(int u, int c){
		long d = 1;
		if(c!=0){
			for(int i=1; i<=c; i++){
				d *= choose(u-i+1,2);
			}
		}
		return d;
	}

	
	/**
	 * The function is to calculate the number of ways that coalescent events on a branch can occur consistently with the gene tree
	 * @param	u	the number of lineages entering the branch
	 * @param	history		coalescent history of the gene tree
	 */
	private long calculateW(TNode node, int c, int[] history){
		long w = 1;
		if(c!=0){
			w = fact(1,c);
			for(int k=0; k<history.length; k++){
				if(history[k]==node.getID()){
					int sum = 0;
					for(int j=0; j<history.length; j++){
						if(j==k || history[j]==-1)continue;
						if((history[j]==node.getID() || S[history[j]][node.getID()]) && R[k][j]){
							sum ++;
						}
					}
					w *= 1.0/(1+sum);
				}
			}
		}
		return w;
	}
	
	
	/**
	 * The function is to calculate the number of ways that coalescent events on a branch can occur consistently with the gene tree divided by c!
	 * It is used for branches under hybridization events.
	 * @param	u	the number of lineages entering the branch
	 * @param	history		coalescent history of the gene tree
	 */
	private double calculateHW(TNode node, int c, int[] history){
		double w = 1;
		if(c!=0){
			for(int k=0; k<history.length; k++){
				if(history[k]==node.getID()){
					int sum = 0;
					for(int j=0; j<history.length; j++){
						if(j==k || history[j]==-1)continue;
						if((history[j]==node.getID() || S[history[j]][node.getID()]) && R[k][j]){
							sum ++;
						}
					}
					w /= 1+sum;
				}
			}
		}
		return w;
	}
	
	
	/**
	 * The function is to do pre-processing for Network2Tree.
	 */
	private void adjustNetwork(Network<Double> net){
		for(NetNode<Double> nnode: net.getNetworkNodes()){
			NetNode<Double> toRemoveNode = nnode.getChildren().iterator().next();
			if(toRemoveNode.isLeaf()){
				continue;
			}
			double distance = toRemoveNode.getParentDistance(nnode);
			for(NetNode<Double> parent: nnode.getParents()){
				nnode.setParentDistance(parent, distance + nnode.getParentDistance(parent));
			}
			int index = 0;
			for(NetNode<Double> child: toRemoveNode.getChildren()){
				distance = child.getParentDistance(toRemoveNode);
			//	nnode.adoptChild(child, distance, toRemoveNode.getGamma(index)*nnode.getGamma(toRemoveNode));
				index++;
			}
			nnode.removeChild(toRemoveNode);
		}
	}
	
	
	/**
	 * The function is to print matrix for debugging
	 */
	private void printMatrix(boolean[][] matrix){
		for(int i=0; i<matrix.length; i++){
			for(int j=0; j<matrix[0].length; j++){
				if(matrix[i][j])
					System.out.print(1 + "\t");
				else
					System.out.print(0 + "\t");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	
	/**
	 * The function is to print a list of coalescent histories for debugging
	 */
	private void printHistories(List<int[]> histories){
		System.out.println("total size:"+histories.size());
		for(int[] his: histories){
			System.out.print("[");
			for(int edge: his){
				System.out.print(edge+" ");
			}
			System.out.println("]");
		}
	}
	
	/**
	 * The function is to print coalescent histories for debugging
	 */
	private void printHistory(int[] history){
		System.out.print("[");
		for(int edge: history){
			System.out.print(edge+" ");
		}
		System.out.println("]");
	}
	
	
	
	/**
	 * The class is to store cluster - edge pair.
	 */
	private class CEPair {
		
		public int cluster_id;
		public int edge_id;
		
		public CEPair(){}
		
		public CEPair(int cluster, int edge) {
			edge_id = edge;
			cluster_id = cluster;
		}
		
		public void set(int cluster, int edge) {
			edge_id = edge;
			cluster_id = cluster;
		}
		
		public int hashCode() {
			return edge_id;
		}
		
		public boolean equals(Object o) {
			if(!(o instanceof CEPair)) {
				return false;
			}
			
			CEPair p2 = (CEPair) o;
			
			return (cluster_id == p2.cluster_id) && (edge_id == p2.edge_id);
		}
		
		public String toString(){
			return "edge:"+edge_id+"/node:"+cluster_id;
		}
	}
	
	
}

