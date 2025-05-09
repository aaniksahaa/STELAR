package phylonet.coalescent;

import java.util.AbstractMap;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.ConcurrentHashMap;

import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;
import phylonet.coalescent.MGDInference_DP.TaxonNameMap;

public class DuplicationWeightCounter {

	HashMap<STBipartition, Integer> weights;
	String[] gtTaxa;
	String[] stTaxa;

	// private List<Set<STBipartition>> X;

	private Map<STBipartition, Integer> geneTreeSTBCount;
	private Map<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer> geneTreeInvalidSTBCont;

	private boolean rooted;

	private TaxonNameMap taxonNameMap;

	private ClusterCollection clusters;

	// Used for dup loss calculations
	List<STITreeCluster> treeAlls = new ArrayList<STITreeCluster>();

	HashMap<STBipartition, Set<STBipartition>> alreadyWeigthProcessed = new HashMap<STBipartition, Set<STBipartition>>();

	public DuplicationWeightCounter(String[] gtTaxa, String[] stTaxa,
			boolean rooted, TaxonNameMap taxonNameMap,
			ClusterCollection clusters2) {
		this.gtTaxa = gtTaxa;
		this.stTaxa = stTaxa;
		this.rooted = rooted;
		this.taxonNameMap = taxonNameMap;
		this.clusters = clusters2;
	}

	private String getSpeciesName(String geneName) {
		String stName = geneName;
		if (taxonNameMap != null) {
			stName = taxonNameMap.getTaxonName(geneName);
		}
		return stName;
	}

	private boolean addToClusters(STITreeCluster c, int size,
			boolean geneTreeCluster) {
		Vertex nv = c.new Vertex();
		return clusters.addCluster(nv, size);
	}

	Integer calculateHomomorphicCost(List<Integer> El, STITreeCluster cluster,
			Vertex smallV, Vertex bigv, List<Tree> trees) {
		Integer e = 0;
		for (int k = 0; k < trees.size(); k++) {
			Tree tr = trees.get(k);
			STITreeCluster treeAll = treeAlls.get(k);
			if (smallV.getCluster().isDisjoint(treeAll)
					|| bigv.getCluster().isDisjoint(treeAll)) {
				continue;
			}
			if (El.get(k) == null) {
				if (taxonNameMap == null) {
					El.set(k, DeepCoalescencesCounter.getClusterCoalNum_rooted(
							tr, cluster));
				} else {
					El.set(k, DeepCoalescencesCounter.getClusterCoalNum_rooted(
							tr, cluster, taxonNameMap));
				}
			}
			e += El.get(k);
		}
		return e;
	}

	/*
	 * Calculates the cost of a cluster based on the standard definition
	 * */
	int calculateDLstdClusterCost(STITreeCluster cluster, List<Tree> trees) {
		/*
		 * if (XLweights.containsKey(cluster)) { return XLweights.get(cluster);
		 * }
		 */
		int weight = 0;
		for (Entry<STBipartition, Integer> entry : geneTreeSTBCount.entrySet()) {
			STBipartition otherSTB = entry.getKey();
			/*if (cluster.containsCluster(otherSTB.c))
				continue;*/
			boolean c1 = cluster.containsCluster(otherSTB.cluster1);
			boolean c2 = cluster.containsCluster(otherSTB.cluster2);
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (Entry<SimpleEntry<STITreeCluster, STITreeCluster>, Integer> entry : geneTreeInvalidSTBCont.entrySet()) {
			SimpleEntry<STITreeCluster, STITreeCluster> otherSTB = entry.getKey();
			boolean c1 = cluster.containsCluster(otherSTB.getKey());
			boolean c2 = cluster.containsCluster(otherSTB.getValue());
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (STITreeCluster treeAll : treeAlls) {
			if (cluster.containsCluster(treeAll)) {
				weight++;
			}
		}
		int ret = weight;
		// XLweights.put(cluster, ret);
		return ret;
	}

	double computeTreeSTBipartitions(MGDInference_DP inference) {

		double unweigthedConstant = 0;
		double weightedConstant = 0;
		int k = inference.trees.size();
		String[] leaves = stTaxa;
		int n = leaves.length;
		boolean duploss = (inference.optimizeDuploss == 3);		

		geneTreeSTBCount = new HashMap<STBipartition, Integer>(k * n);
		geneTreeInvalidSTBCont = new HashMap<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer>();
		// geneTreeRootSTBs = new HashMap<STBipartition, Integer>(k*n);
		// needed for fast version
		// clusterToSTBs = new HashMap<STITreeCluster, Set<STBipartition>>(k*n);

		STITreeCluster all = new STITreeCluster(stTaxa);
		String as[];
		int j = (as = stTaxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i];
			all.addLeaf(t);
		}
		addToClusters(all, leaves.length, false);

		for (int t = 0; t < inference.trees.size(); t++) {
			Tree tr = inference.trees.get(t);

			STITreeCluster allInducedByGT = new STITreeCluster(stTaxa);

			String[] gtLeaves = tr.getLeaves();
			for (int i = 0; i < gtLeaves.length; i++) {
				String l = gtLeaves[i];
				allInducedByGT.addLeaf(getSpeciesName(l));
			}
			treeAlls.add(allInducedByGT);
			int allInducedByGTSize = allInducedByGT.getClusterSize();

			weightedConstant += duploss ? 
					2 * (allInducedByGTSize - 1) : 0;
					
			unweigthedConstant += (tr.getLeafCount() - 1);

			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
			// Map<TNode,STITreeCluster> nodeToGTCluster = new HashMap<TNode,
			// STITreeCluster>(n);

			for (TNode node : tr.postTraverse()) {				
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					String nodeName = node.getName();

					// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
					// gtCluster.addLeaf(nodeName);

					nodeName = getSpeciesName(nodeName);
					STITreeCluster cluster = new STITreeCluster(leaves);
					cluster.addLeaf(nodeName);

					addToClusters(cluster, 1, true);

					nodeToSTCluster.put(node, cluster);

					// nodeToGTCluster.put(node, gtCluster);

					if (!rooted) {
						throw new RuntimeException("Unrooted not implemented.");
						// tryAddingSTB(cluster, treeComplementary(null
						// /*gtCluster*/,leaves), null, node, true);
					}
				} else {
					int childCount = node.getChildCount();
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(leaves.length);
					// BitSet gbs = new BitSet(leaves.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
						// gbs.or(nodeToGTCluster.get(child).getBitSet());
					}

					// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
					// gtCluster.setCluster(gbs);

					STITreeCluster cluster = new STITreeCluster(leaves);
					cluster.setCluster((BitSet) bs.clone());

					int size = cluster.getClusterSize();

					addToClusters(cluster, size, true);
					nodeToSTCluster.put(node, cluster);
					// nodeToGTCluster.put(node, gtCluster);

					if (rooted) {

						if (index > 2) {
							throw new RuntimeException(
									"None bifurcating tree: " + tr + "\n"
											+ node);
						}

						STITreeCluster l_cluster = childbslist[0];
						STITreeCluster r_cluster = childbslist[1];

						tryAddingSTB(l_cluster, r_cluster, cluster, node, true);

					} else {
						throw new RuntimeException("Unrooted not implemented.");
						/*
						 * if (childCount == 2) { STITreeCluster l_cluster =
						 * childbslist[0];
						 * 
						 * STITreeCluster r_cluster = childbslist[1];
						 * 
						 * STITreeCluster allMinuslAndr_cluster =
						 * treeComplementary(null this should be
						 * gtCluster?,leaves);
						 * 
						 * STITreeCluster lAndr_cluster = cluster;
						 * 
						 * if (allMinuslAndr_cluster.getClusterSize() != 0) { //
						 * add Vertex STBs tryAddingSTB(l_cluster, r_cluster,
						 * cluster, node, true); tryAddingSTB( r_cluster,
						 * allMinuslAndr_cluster, null, node, true);
						 * tryAddingSTB(l_cluster, allMinuslAndr_cluster, null,
						 * node, true);
						 * 
						 * // Add the Edge STB tryAddingSTB(lAndr_cluster,
						 * allMinuslAndr_cluster, null, node, true); }
						 * 
						 * } else if (childCount == 3 && node.isRoot()) {
						 * STITreeCluster l_cluster = childbslist[0];
						 * 
						 * STITreeCluster m_cluster = childbslist[1];
						 * 
						 * STITreeCluster r_cluster = childbslist[2];
						 * 
						 * tryAddingSTB(l_cluster, r_cluster, null, node, true);
						 * tryAddingSTB(r_cluster, m_cluster, null, node, true);
						 * tryAddingSTB(l_cluster, m_cluster, null, node, true);
						 * } else { throw new
						 * RuntimeException("None bifurcating tree: "+ tr+ "\n"
						 * + node); }
						 */}
				}
			}

		}

		int s = 0;
		for (Integer c : geneTreeSTBCount.values()) {
			s += c;
		}
		System.err.println("STBs in gene trees (count): "
				+ geneTreeSTBCount.size());
		System.err.println("STBs in gene trees (sum): " + s);

		s = clusters.getClusterCount();

		System.err.println("Number of Clusters: " + s);

		weights = new HashMap<STBipartition, Integer>(
				geneTreeSTBCount.size() * 2);
		// System.err.println("sigma n is "+sigmaN);

		if (inference.DLbdWeigth == -1) {			
			inference.DLbdWeigth = (weightedConstant + 2*k + 0.0D) / 2*(k*n);
			System.out.println("Estimated bd weight = " + inference.DLbdWeigth);
		}
			
		return (unweigthedConstant + (1 - inference.DLbdWeigth) * weightedConstant);
	}

	void addAllPossibleSubClusters(STITreeCluster cluster, int size) {
		BitSet bs = (BitSet) cluster.getBitSet().clone();
		bs.clear(0, size);
		while (true) {
			int tsb = bs.nextClearBit(0);
			if (tsb >= size) {
				break;
			}
			bs.set(tsb);
			bs.clear(0, tsb);
			STITreeCluster c = new STITreeCluster(cluster.getTaxa());
			c.setCluster((BitSet) bs.clone());
			addToClusters(c, c.getClusterSize(), false);
		}
		System.err
				.println("Number of Clusters After Adding All possible clusters: "
						+ clusters.getClusterCount());
	}

	void addExtraBipartitionsByInput(ClusterCollection extraClusters,
			List<Tree> trees, boolean extraTreeRooted) {

		String[] leaves = stTaxa;
		int n = leaves.length;

		// STITreeCluster all = extraClusters.getTopVertex().getCluster();

		for (Tree tr : trees) {

			/*
			 * String[] treeLeaves = tr.getLeaves(); STITreeCluster treeAll =
			 * new STITreeCluster(treeLeaves); for (int i = 0; i <
			 * treeLeaves.length; i++) { String l = treeLeaves[i];
			 * treeAll.addLeaf(l); }
			 */
			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(
					n);

			for (Iterator<TNode> nodeIt = tr.postTraverse().iterator(); nodeIt
					.hasNext();) {
				TNode node = nodeIt.next();
				if (node.isLeaf()) {
					String treeName = node.getName();
					String nodeName = getSpeciesName(treeName);

					STITreeCluster tb = new STITreeCluster(leaves);
					tb.addLeaf(nodeName);

					nodeToSTCluster.put(node, tb);

					if (!extraTreeRooted) {
						// TODO: fix the following (first null)
						throw new RuntimeException("Unrooted not implemented.");
						// tryAddingSTB(tb, tb.complementaryCluster(),
						// all,node,false);
					}
				} else {
					int childCount = node.getChildCount();
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(leaves.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster(leaves);
					cluster.setCluster((BitSet) bs.clone());

					addToClusters(cluster, cluster.getClusterSize(), false);
					nodeToSTCluster.put(node, cluster);

					if (extraTreeRooted) {

						if (index > 2) {
							throw new RuntimeException(
									"None bifurcating tree: " + tr + "\n"
											+ node);
						}

						STITreeCluster l_cluster = childbslist[0];

						STITreeCluster r_cluster = childbslist[1];

						tryAddingSTB(l_cluster, r_cluster, cluster, node, false);
					} else {
						throw new RuntimeException("Unrooted not implemented.");
						/*
						 * if (childCount == 2) { STITreeCluster l_cluster =
						 * childbslist[0];
						 * 
						 * STITreeCluster r_cluster = childbslist[1]; // Fix the
						 * following (first null) STITreeCluster
						 * allMinuslAndr_cluster = treeComplementary(null,
						 * leaves);
						 * 
						 * STITreeCluster lAndr_cluster = cluster;
						 * 
						 * // add Vertex STBs tryAddingSTB( l_cluster,
						 * r_cluster, cluster,node,false); if
						 * (allMinuslAndr_cluster.getClusterSize() != 0) {
						 * tryAddingSTB( r_cluster, allMinuslAndr_cluster,
						 * null,node,false); tryAddingSTB( l_cluster,
						 * allMinuslAndr_cluster, null,node,false);
						 * tryAddingSTB( lAndr_cluster, allMinuslAndr_cluster,
						 * all,node,false); }
						 * 
						 * } else if (childCount == 3 && node.isRoot()) {
						 * STITreeCluster l_cluster = childbslist[0];
						 * 
						 * STITreeCluster m_cluster = childbslist[1];
						 * 
						 * STITreeCluster r_cluster = childbslist[2];
						 * 
						 * tryAddingSTB( l_cluster, r_cluster, null,node,false);
						 * tryAddingSTB( r_cluster, m_cluster, null,node,false);
						 * tryAddingSTB( l_cluster, m_cluster, null,node,false);
						 * } else { throw new
						 * RuntimeException("None bifurcating tree: "+ tr+ "\n"
						 * + node); }
						 */
					}
				}
			}

		}
		int s = extraClusters.getClusterCount();
		/*
		 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
		 */
		System.err
				.println("Number of Clusters After additions from extra Trees: "
						+ s);
	}

	private void tryAddingSTB(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster cluster, TNode node,
			boolean fromGeneTrees) {
		// System.err.println("before adding: " + STBCountInGeneTrees);
		// System.err.println("Trying: " + l_cluster + "|" + r_cluster);
		int size = cluster.getClusterSize();
		if (l_cluster.isDisjoint(r_cluster)) {

			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);
			((STINode) node).setData(stb);
			if (fromGeneTrees) {
				clusters.addGeneTreeSTB(stb, size);
				// gtNodeToSTBs.put(node,stb);
				// addSTBToX(stb,size);
				// System.out.println(stb + " hashes to " + stb.hashCode());
				// if (! hash.containsKey(stb.hashCode()))
				// hash.put(stb.hashCode(), new HashSet<STBipartition>());
				// hash.get(stb.hashCode()).add(stb);
				geneTreeSTBCount.put(
						stb,
						geneTreeSTBCount.containsKey(stb) ? geneTreeSTBCount
								.get(stb) + 1 : 1);
			}

			/*
			 * if (size == allInducedByGTSize){ if (!
			 * geneTreeRootSTBs.containsKey(stb)) { geneTreeRootSTBs.put(stb,
			 * 1); } else { geneTreeRootSTBs.put(stb,
			 * geneTreeRootSTBs.get(stb)+1); } }
			 */
		} else {
			AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster> stb = 
					l_cluster.getBitSet().cardinality() > r_cluster.getBitSet().cardinality() ? 
					new AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>(l_cluster, r_cluster) : 
					new AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>(r_cluster, l_cluster);
			geneTreeInvalidSTBCont.put(stb,	geneTreeInvalidSTBCont.containsKey(stb) ? 
					geneTreeInvalidSTBCont.get(stb) + 1 : 
					1);
			// System.err.println("Adding only to extra");
			// This case could happen for multiple-copy
			BitSet and = (BitSet) l_cluster.getBitSet().clone();
			and.and(r_cluster.getBitSet());

			BitSet l_Minus_r = (BitSet) and.clone();
			l_Minus_r.xor(l_cluster.getBitSet());
			STITreeCluster lmr = new STITreeCluster(stTaxa);
			lmr.setCluster(l_Minus_r);

			BitSet r_Minus_l = (BitSet) and.clone();
			r_Minus_l.xor(r_cluster.getBitSet());
			STITreeCluster rml = new STITreeCluster(stTaxa);
			rml.setCluster(r_Minus_l);

			if (!rml.getBitSet().isEmpty()) {
				addToClusters(rml, rml.getClusterSize(), false);
				// addSTBToX( new STBipartition(l_cluster, rml, cluster),size);
			}
			if (!lmr.getBitSet().isEmpty()) {
				addToClusters(lmr, lmr.getClusterSize(), false);
				// addSTBToX(new STBipartition(lmr, r_cluster, cluster), size);
			}
		}
	}

	public Integer getCalculatedBiPartitionDPWeight(STBipartition bi) {
		if (!weights.containsKey(bi)) {
			// weights.put(bi,calculateMissingWeight(bi));
			return null;
		}
		return weights.get(bi);
	}

	// static public int cnt = 0;

	void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {

		if (rooted && taxonNameMap == null && stTaxa.length > trees.size()) {
			//calculateWeightsByLCA(trees, trees);
			if (extraTrees != null) {
				//calculateWeightsByLCA(extraTrees, trees);
			}
		}

	}

	void calculateWeightsByLCA(List<Tree> stTrees, List<Tree> gtTrees) {

		for (Tree stTree : stTrees) {
			SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(stTree);
			for (Tree gtTree : gtTrees) {
				Stack<TNode> stack = new Stack<TNode>();
				for (TNode gtNode : gtTree.postTraverse()) {
					if (gtNode.isLeaf()) {
						stack.push(stTree.getNode(gtNode.getName()));
					} else {
						TNode rightLCA = stack.pop();
						TNode leftLCA = stack.pop();
						// If gene trees are incomplete, we can have this case
						if (rightLCA == null || leftLCA == null) {
							stack.push(null);
							continue;
						}
						TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
						stack.push(lca);
						if (lca != leftLCA && lca != rightLCA) {
							// LCA in stTree dominates gtNode in gene tree
							// gtTree
							STBipartition stSTB = (STBipartition) ((STINode) lca)
									.getData();
							STBipartition gtSTB = (STBipartition) ((STINode) gtNode)
									.getData();
							Set<STBipartition> alreadyProcessedSTBs = alreadyWeigthProcessed
									.get(gtSTB);

							if (alreadyProcessedSTBs == null) {
								alreadyProcessedSTBs = new HashSet<STBipartition>(
										gtTrees.size() / 4);
								alreadyWeigthProcessed.put(gtSTB,
										alreadyProcessedSTBs);
							}

							if (alreadyProcessedSTBs.contains(stSTB)) {
								continue;
							}

							weights.put(
									stSTB,
									(weights.containsKey(stSTB) ? weights
											.get(stSTB) : 0)
											+ geneTreeSTBCount.get(gtSTB));
							alreadyProcessedSTBs.add(stSTB);
						}
					}
				}
			}
		}
	}

	class CalculateWeightTask {
		
		// Inner class to store matching counts for X and Y taxa
	    private static class MatchingCounts {
	        int xCount;
	        int yCount;
	        
	        MatchingCounts(int x, int y) {
	            this.xCount = x;
	            this.yCount = y;
	        }
	        
	        void add(MatchingCounts other) {
	            this.xCount += other.xCount;
	            this.yCount += other.yCount;
	        }
	    }

		/**
		 * 
		 */
		private static final long serialVersionUID = -2614161117603289345L;
		private STBipartition stb;
		private ClusterCollection containedClusterCollection;
		private List<Tree> trees;
		
		
		// Utility class to store BitSet pairs with proper equals & hashCode
	    private static class BitSetPair {
	        private final BitSet first;
	        private final BitSet second;
	        private final int hash;

	        BitSetPair(BitSet a, BitSet b) {
	            // Ensure consistent order
	            if (a.hashCode() <= b.hashCode()) {
	                this.first = (BitSet) a.clone();
	                this.second = (BitSet) b.clone();
	            } else {
	                this.first = (BitSet) b.clone();
	                this.second = (BitSet) a.clone();
	            }
	            this.hash = first.hashCode() ^ second.hashCode();
	        }

	        @Override
	        public boolean equals(Object obj) {
	            if (this == obj) return true;
	            if (!(obj instanceof BitSetPair)) return false;
	            BitSetPair other = (BitSetPair) obj;
	            return first.equals(other.first) && second.equals(other.second);
	        }

	        @Override
	        public int hashCode() {
	            return hash;
	        }
	    }
		
		// Cache to store computed values for (X, Y)
	    private static final Map<BitSetPair, Integer> weightCache = new ConcurrentHashMap<>();
		

		public CalculateWeightTask(STBipartition stb,
				ClusterCollection collection, List<Tree>trees) {
			this.stb = stb;
			this.containedClusterCollection = collection;
			this.trees = trees;
		}
		
		int calculateMissingWeight() {
			//System.err.print("Calculating weight for: " + biggerSTB);
//			int weight = 0;
			
			BitSet X =  stb.cluster1.getBitSet() ;
			BitSet Y =  stb.cluster2.getBitSet() ;
			
			// Normalize and create a key
	        BitSetPair key = new BitSetPair(X, Y);
	        
	        boolean isWeightCaching = true;
	        isWeightCaching = false;

	        // Check cache
	        // this caching is, in turn, increasing runtime, not much helping
	        
	        if(isWeightCaching) {
		        return weightCache.computeIfAbsent(key, k -> {
		            // Perform computation if not cached
		            return computeWeight();
		        });
	        } else {
	        	return computeWeight();
	        }
	        
		}

		int computeWeight() {
			//System.err.print("Calculating weight for: " + biggerSTB);
			int weight = 0;
			
			BitSet X =  stb.cluster1.getBitSet() ;
			BitSet Y =  stb.cluster2.getBitSet() ;
			//System.out.println("Calculating missing Weight");
			//System.out.println(stb.toString());
			
//			System.out.println("\nTaxa:\n");
//			for(String t: stb.cluster1.getTaxa()) {
//				System.out.println("Taxon " + t + " ");
//			}
//			System.out.println("\n");
			
			// traversing all gene trees
//			System.out.println("total trees = " + trees.size());
			
			
			
			int w1 = 0;
			int w2 = 0;
			
			// isNew variable controls whether the improved or the previous algorithm is being run
			boolean isNew = false;
//			isNew = true;
			
			// O(nk) implementation instead of O(n^2k)
			if(isNew) {
				int cntt = 0;
				String[] Xtaxa = stb.cluster1.getTaxa();
				String[] Ytaxa = stb.cluster2.getTaxa();
				
				Set<String> xTaxaSet = new HashSet<>();
		        Set<String> yTaxaSet = new HashSet<>();
		        
		        for(int i=0; i<Xtaxa.length; i++) {
		        	if(X.get(i)) {
		        		xTaxaSet.add(Xtaxa[i]);
		        	}
		        }
		        
		        for(int i=0; i<Ytaxa.length; i++) {
		        	if(Y.get(i)) {
		        		yTaxaSet.add(Ytaxa[i]);
		        	}
		        }
		        
		        
		        // /* prev imp
				
		        // looping over the gene trees
				for (int t = 0; t < trees.size(); t++) {
					Tree tr = trees.get(t);
					
					// for each node, we need to keep two numbers
					// matching leaves count (in its subtree) with X and Y
					// Map<TNode, MatchingCounts> nodeCounts = new HashMap<>();
					
					// this stack actually does the DFS
					// here what we do is, say we have a node A with children L and R
	                // and let, for node N, N.x and N.y denote the match counts with X and Y respectively
	                // MatchingCounts class serves this purpose of a pair of numbers
					Stack<MatchingCounts> stack = new Stack<>();
					
					// we postTraverse, since it corresponds to DFS
					for (TNode node : tr.postTraverse()) {
						// System.err.println("Node is:" + node);
						if (node.isLeaf()) {
							// from nodeName determine its matching
							// since this is just one node, the count should be 0 or 1
							String nodeName = node.getName();
							
							
							// this id is not working, need to understand more
//							int nodeID = node.getID() - 1;
//							int xMatch = X.get(nodeID) ? 1 : 0;
//			                int yMatch = Y.get(nodeID) ? 1 : 0;
							
							
							// Check if the leaf matches with X or Y taxa
			                
			                int xMatch = xTaxaSet.contains(nodeName) ? 1 : 0;
			                int yMatch = yTaxaSet.contains(nodeName) ? 1 : 0;
			                
			                MatchingCounts mc = new MatchingCounts(xMatch, yMatch);
			                
			                //nodeCounts.put(node, mc);
			                
			                stack.push(mc);
						} else {
							cntt++;
							// For internal nodes, accumulate counts from children
			                MatchingCounts currentCounts = new MatchingCounts(0, 0);
			                
//			                MatchingCounts L = null, R = null;
//			                int i = 0;
			                
			                // here what we do is, say we have a node A with children L and R
			                // and let, for node N, N.x and N.y denote the match counts with X and Y respectively
			                // thus, we can say, A.x = L.x + R.x and A.y = L.y + R.y
			                // since we are using a stack, we may just add the counts after popping the children off the stack
			                MatchingCounts R = stack.pop();
		                	MatchingCounts L = stack.pop();
		                	currentCounts.add(L);
		                	currentCounts.add(R);
		                	stack.push(currentCounts);
		                	
//			                for (TNode child : node.getChildren()) {
//			                	
//			                    MatchingCounts childCounts = nodeCounts.get(child);
//			                    if(i == 0) {
//			                    	L = childCounts;
//			                    } else {
//			                    	R = childCounts;
//			                    }
//			                    if (childCounts != null) {
//			                        currentCounts.add(childCounts);
//			                    }
//			                    i++;
//			                }
			                
			       
			                //nodeCounts.put(node, currentCounts);
			                
		                	// here, after we have the X and Y matches, we now add the triplet scores up the node
			                int temp = 0, res;
			                
			                int c1 = L.xCount;
			        		int c2 = R.yCount;
	
			        		res = c1*c2*(c1+c2-2);
			        		res /=2;
			        		
			        		temp += res;
			        		
			        		c1 = L.yCount;
			        		c2 = R.xCount;
	
			        		res = c1*c2*(c1+c2-2);
			        		res /=2;
			        		
			        		temp += res;
				
							//System.out.println(temp);
			        		
			        		// finally we accumate it to the final total weight count means the total triplet score
			        		weight += temp;
			        		
						}
					}
				}
				w1 = weight;
				//System.out.println("Total "+cntt+" STBs");
				 
		        
		        
		        
		        /* gpt improve imp
		        
		        // looping over the gene trees
		        for (int t = 0; t < trees.size(); t++) {
		            Tree tr = trees.get(t);
		            // Use an ArrayDeque as a stack for better performance
		            Deque<int[]> stack = new ArrayDeque<>();
		            
		            // Perform post-order traversal
		            for (TNode node : tr.postTraverse()) {
		                if (node.isLeaf()) {
//		                    int nodeID = node.getID();
//		                    int xMatch = X.get(nodeID) ? 1 : 0;
//		                    int yMatch = Y.get(nodeID) ? 1 : 0;
		                    
		                    
		                    String nodeName = node.getName();
		                    int xMatch = xTaxaSet.contains(nodeName) ? 1 : 0;
			                int yMatch = yTaxaSet.contains(nodeName) ? 1 : 0;
		                    
		                    
		                    // Create an array to hold the match counts for this leaf node
		                    int[] counts = new int[2];
		                    counts[0] = xMatch;  // x count
		                    counts[1] = yMatch;  // y count
		                    stack.push(counts);
		                } else {
		                    // Pop counts for the two children
		                    int[] right = stack.pop();
		                    int[] left = stack.pop();
		                    
		                    // Aggregate counts for the current internal node
		                    int[] current = new int[2];
		                    current[0] = left[0] + right[0]; // Sum of x counts
		                    current[1] = left[1] + right[1]; // Sum of y counts
		                    stack.push(current);
		                    
		                    // Calculate the triplet score for the internal node
		                    int temp = 0;
		                    // Case 1: left.x * right.y triplets
		                    int c1 = left[0];
		                    int c2 = right[1];
		                    int res = c1 * c2 * (c1 + c2 - 2) / 2;
		                    temp += res;
		                    
		                    // Case 2: left.y * right.x triplets
		                    c1 = left[1];
		                    c2 = right[0];
		                    res = c1 * c2 * (c1 + c2 - 2) / 2;
		                    temp += res;
		                    
		                    // Accumulate the computed triplet score
		                    weight += temp;
		                }
		            }
		        }
		        w1 = weight;
		        
		        *
		        */

			}
			else {   
				
				// the following is the previous implementation
				// Note that, this does one clever thing, rather than looping over all trees
				// it has a pre-determined aggregated list of all the unique STBs in the gene trees
				// basically it has a list of unique STBs and also a count map
				// so, it calculates for each of them separately, and multiplies by count of appearance of that STB, then adds to the total
				// since, practically there are many overlaps of STBs, this improves the performance a lot against the theoretical bound of O(n^2k)
				
				weight = 0;
				
				int cntt = 0;
				for (STBipartition smallerSTB : clusters.getContainedGeneTreeSTBs()) {
					
						int temp = 0;
						// possible place of implementation.
					//	System.out.println("Loop "+ smallerSTB.toString() +"count =  "+ geneTreeSTBCount.get(smallerSTB)); 
						BitSet A = smallerSTB.cluster1.getBitSet();
						BitSet B = smallerSTB.cluster2.getBitSet();
						
						BitSet X1 = and(X,A);
						BitSet Y1 = and(Y,B);
						temp = apply(X1, Y1);
	
						BitSet X2 = and(X,B);
						BitSet Y2 = and(Y,A);
						temp += apply(X2, Y2);
									
						//System.out.println(geneTreeSTBCount.get(smallerSTB));
						temp *= geneTreeSTBCount.get(smallerSTB);
						weight += temp;
						//System.out.println(smallerSTB.toString() + " :: "+ temp);
						
						cntt += geneTreeSTBCount.get(smallerSTB);
				}
				//System.out.println("Total "+cntt+" STBs");
				
				w2 = weight;
			}
			
//			if(w1 != w2) {
//				System.out.println("\nNot Equal!! w1 = " + w1 + ", w2 = " + w2);
//			}
						
		//	System.out.println("STB score of + " + stb.toString() + " is =  "+weight);
			// System.err.print(" ... " + weight);
			
			
			
			// this is the place where changes are meant to be made in case of unrooted STELAR
			
			if (!rooted) {
				throw new RuntimeException("Unrooted not implemented.");
				/*
				 * for (STBipartition rootSTB : geneTreeRootSTBs.keySet()) { int
				 * c = geneTreeRootSTBs.get(rootSTB); STBipartition inducedSTB =
				 * biggerSTB.getInducedSTB(rootSTB.c); if
				 * (inducedSTB.equals(rootSTB)){ weight -= 2 * c;
				 * //System.err.print(" .. (" + rootSTB +" )" +c+" "+ weight);
				 * if (inducedSTB.cluster1.getClusterSize() != 1 &&
				 * inducedSTB.cluster2.getClusterSize() != 1) { weight -= 2 * c;
				 * //System.err.print(" . " + weight); } } }
				 */
			}
			weights.put(stb, weight);
			// System.err.println("Weight of " + biggerSTB + " is " + weight);
			return weight;
		}

		protected Integer compute() {
			return calculateMissingWeight();
		}

	}
	int apply(BitSet x, BitSet y) {
		int res = 0;
		int c1 = x.cardinality();
		int c2 = y.cardinality();

		res = c1*(c1-1)*c2 + c2*(c2-1)*c1;
		res /=2;
		return res;
	}
	BitSet and(BitSet x, BitSet y) {
		BitSet _x = (BitSet) x.clone();
		BitSet _y = (BitSet) y.clone();
		_x.and(_y);

		return _x;
	}

	/*
	 * public void addGoodSTB (STBipartition good, int size) {
	 * goodSTBs.get(size).add(good); }
	 */
	/*
	 * public Set<STBipartition> getClusterBiPartitions(STITreeCluster cluster)
	 * {
	 * 
	 * return clusterToSTBs.get(cluster); }
	 */

	/*
	 * public boolean addCompleteryVertx(Vertex x, STITreeCluster refCluster) {
	 * STITreeCluster c = x._cluster; Vertex reverse = new Vertex();
	 * reverse._cluster = new STITreeCluster(refCluster);
	 * reverse._cluster.getCluster().xor(c.getCluster()); int size =
	 * reverse._cluster.getClusterSize(); if
	 * (!clusters.get(size).contains(reverse)){ clusters.get(size).add(reverse);
	 * return true; //System.err.println("Clusters: "+clusters); } return false;
	 * }
	 */
	/*
	 * private void addSTBToX(STBipartition stb, int size) {
	 * //System.err.println("Adding to X: "+stb+" "+stb.c
	 * +" "+clusterToSTBs.containsKey(stb.c)); if
	 * (clusterToSTBs.containsKey(stb.c) &&
	 * clusterToSTBs.get(stb.c).contains(stb)){ return; } //int size =
	 * stb.c.getClusterSize(); // TODO: following line is algorithmically
	 * harmless, // but inefficient. is it necessary?
	 * //geneTreeSTBCount.put(stb, 0); addToClusters(stb.c, size, false); //
	 * Following needed for Fast //Set<STBipartition> stbs =
	 * clusterToSTBs.get(c); //stbs = (stbs== null)? new
	 * HashSet<STBipartition>() : stbs; //stbs.add(stb); //clusterToSTBs.put(c,
	 * stbs); //System.err.println("X updated: "+STBCountInGeneTrees);
	 * //System.err.println("X updated: "+clusterToSTBs); }
	 */

	/*
	 * private STITreeCluster treeComplementary(STITreeCluster gtCluster,
	 * String[] leaves){ //System.err.print("Tree complementary of "+gtCluster);
	 * STITreeCluster newGTCluster = gtCluster.complementaryCluster();
	 * //System.err.println(" is: "+newGTCluster.getCluster()); STITreeCluster
	 * newSTCluster = new STITreeCluster(leaves); for (String s :
	 * newGTCluster.getClusterLeaves()) {
	 * newSTCluster.addLeaf(getSpeciesName(s)); }
	 * //System.err.println("Tree complementary of "
	 * +gtCluster+" is: "+newSTCluster); return newSTCluster; }
	 */

	/*
	 * private STITreeCluster treeComplementary(List<String> treeNames, Cluster
	 * c , TaxonNameMap taxonMap){ HashSet<String> set = new HashSet<String> ();
	 * set.add(cluster); return treeComplementary(treeNames, set, taxonMap); }
	 */

	/*
	 * void addExtraBipartitionsByHeuristics(ClusterCollection clusters2) {
	 * //goodSTBs = X; //if (true) return; int added = 0; for (int i=1;
	 * i<goodSTBs.size(); i++) { Set<STBipartition> curr_set = goodSTBs.get(i);
	 * for (STBipartition stb1:curr_set) { //if (Math.random() < 0.70) continue;
	 * for (int j=i; j<goodSTBs.size(); j++) { Set<STBipartition> other_set =
	 * goodSTBs.get(j); //if (Math.random() < 0.70) continue; for (STBipartition
	 * stb2:other_set) { //System.out.println(stb1 +" **AND** " + stb2); if
	 * (stb1.cluster1.getClusterSize() < 3 || stb1.cluster2.getClusterSize() < 3
	 * || stb2.cluster1.getClusterSize() < 3 || stb2.cluster2.getClusterSize() <
	 * 3) { if (tryToAdd(stb1,stb2,bipToAddToX) != null) added++; } }
	 * System.err.println(bipToAddToX.size() + " " + i); } } }
	 * 
	 * for (STBipartition stb: bipToAddToX) { //System.err.println( "Adding: " +
	 * stb); addSTBToX(clusters, stb); } System.out.println("\n\nAdded " +
	 * added+ " bipartitions:\n");
	 * 
	 * int s = 0; for (Integer c: clusters.keySet()){ s +=
	 * clusters.get(c).size(); }
	 * System.out.println("Number of Clusters After Addition: " +s);
	 * 
	 * }
	 * 
	 * 
	 * private STBipartition tryAddingExtraSTB_AndreRule(STBipartition stb1,
	 * STBipartition stb2, Set<STBipartition> bipToAddToX) { if
	 * (stb1.equals(stb2)) return null; if ( stb1.isDominatedBy(stb2) ||
	 * stb2.isDominatedBy(stb1) ) return null;
	 * 
	 * if ( stb1.c.isDisjoint(stb2.c) ) return null;
	 * 
	 * if ( stb1.cluster1.isDisjoint(stb2.cluster2) &&
	 * stb1.cluster2.isDisjoint(stb2.cluster1)) { STITreeCluster cl1 = new
	 * STITreeCluster(stb1.cluster1); cl1 = cl1.merge(stb2.cluster1);
	 * STITreeCluster cl2 = new STITreeCluster(stb1.cluster2); cl2 =
	 * cl2.merge(stb2.cluster2); STITreeCluster cl = new STITreeCluster(stb1.c);
	 * cl = cl.merge(stb2.c); STBipartition r = new STBipartition(cl1,cl2,cl);
	 * bipToAddToX.add(r); return r; } else if (
	 * stb1.cluster1.isDisjoint(stb2.cluster1) &&
	 * stb1.cluster2.isDisjoint(stb2.cluster2) ) { STITreeCluster cl1 = new
	 * STITreeCluster(stb1.cluster1); cl1 = cl1.merge(stb2.cluster2);
	 * STITreeCluster cl2 = new STITreeCluster(stb1.cluster2); cl2 =
	 * cl2.merge(stb2.cluster1); STITreeCluster cl = new STITreeCluster(stb1.c);
	 * cl = cl.merge(stb2.c); STBipartition r = new STBipartition(cl1,cl2,cl);
	 * bipToAddToX.add(r); return r; } return null; }
	 */
}
