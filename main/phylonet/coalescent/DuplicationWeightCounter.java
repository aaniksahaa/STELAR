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

	HashMap<STBipartition, Double> weights;
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

		weights = new HashMap<STBipartition, Double>(
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

	public Double getCalculatedBiPartitionDPWeight(STBipartition bi) {
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
	    private static class DPValuesUnweighted {
	        double PX;
	        double PY;
	        
	        DPValuesUnweighted(double PX, double PY) {
	        	this.PX = PX;
	            this.PY = PY;
	        }
	    }
	    
	    DPValuesUnweighted getParentDPValues(DPValuesUnweighted L, DPValuesUnweighted R) {
	        double PX_A = L.PX + R.PX;
	        double PY_A = L.PY + R.PY;
	        
	        // Return new DPValues for the parent node A (parentDistance set to 0 since it's not specified)
	        return new DPValuesUnweighted(PX_A, PY_A);
	    }
	    
	    private static class DPValuesWeighted1 {
	        double PX;
	        double PY;
	        double QX;
	        double QY;
	        double parentDistance;
	        
	        DPValuesWeighted1(double PX, double PY, double QX, double QY, double parentDistance) {
	            this.PX = PX;
	            this.PY = PY;
	            this.QX = QX;
	            this.QY = QY;
	            this.parentDistance = parentDistance;
	        }
	    }
	    
	    
	    
	    DPValuesWeighted1 getParentDPValues(TNode parentNode, DPValuesWeighted1 L, DPValuesWeighted1 R) {
	        // Compute the exponential decay factors
	        double expWAB = Math.exp(-L.parentDistance); // e^(-w(A,B))
	        double expWAC = Math.exp(-R.parentDistance); // e^(-w(A,C))
	        
	        // Compute PX(A) = e^(-w(A,B)) * PX(B) + e^(-w(A,C)) * PX(C)
	        double PX_A = expWAB * L.PX + expWAC * R.PX;
	        
	        // Compute QX(A) = e^(-w(A,B)) * QX(B) + e^(-w(A,C)) * QX(C) + e^(-w(A,B)) * PX(B) * e^(-w(A,C)) * PX(C)
	        double QX_A = expWAB * L.QX + expWAC * R.QX + (expWAB * L.PX) * (expWAC * R.PX);
	        
	        // Compute PY(A) = e^(-w(A,B)) * PY(B) + e^(-w(A,C)) * PY(C)
	        double PY_A = expWAB * L.PY + expWAC * R.PY;
	        
	        // Compute QY(A) = e^(-w(A,B)) * QY(B) + e^(-w(A,C)) * QY(C) + e^(-w(A,B)) * PY(B) * e^(-w(A,C)) * PY(C)
	        double QY_A = expWAB * L.QY + expWAC * R.QY + (expWAB * L.PY) * (expWAC * R.PY);
	        
	        // Return new DPValues for the parent node A (parentDistance set to 0 since it's not specified)
	        return new DPValuesWeighted1(PX_A, PY_A, QX_A, QY_A, parentNode.getParentDistance());
	    }
	    
	    
	    
	    
	    private static class DPValuesWeighted2 {
	        double PX;
	        double PY;
	        double QX;
	        double QY;
	        double RX; // New field for R_X
	        double RY; // New field for R_Y
	        double parentDistance;
	        
	        DPValuesWeighted2(double PX, double PY, double QX, double QY, double RX, double RY, double parentDistance) {
	            this.PX = PX;
	            this.PY = PY;
	            this.QX = QX;
	            this.QY = QY;
	            this.RX = RX;
	            this.RY = RY;
	            this.parentDistance = parentDistance;
	        }
	    }
	    
	    DPValuesWeighted2 getParentDPValues(TNode parentNode, DPValuesWeighted2 L, DPValuesWeighted2 R) {
	        // Compute the exponential decay factors
	        double expWAB = Math.exp(-L.parentDistance); // e^(-w(A,B))
	        double expWAC = Math.exp(-R.parentDistance); // e^(-w(A,C))
	        
	        // Compute PX(A) = e^(-w(A,B)) * PX(B) + e^(-w(A,C)) * PX(C)
	        double PX_A = expWAB * L.PX + expWAC * R.PX;
	        
	        // Compute QX(A) = QX(B) + QX(C) + e^(-w(A,B)) * PX(B) * e^(-w(A,C)) * PX(C)
	        double QX_A = L.QX + R.QX + (expWAB * L.PX) * (expWAC * R.PX);
	        
	        // Compute RX(A) = RX(B) + RX(C)
	        double RX_A = L.RX + R.RX;
	        
	        // Compute PY(A) = e^(-w(A,B)) * PY(B) + e^(-w(A,C)) * PY(C)
	        double PY_A = expWAB * L.PY + expWAC * R.PY;
	        
	        // Compute QY(A) = QY(B) + QY(C) + e^(-w(A,B)) * PY(B) * e^(-w(A,C)) * PY(C)
	        double QY_A = L.QY + R.QY + (expWAB * L.PY) * (expWAC * R.PY);
	        
	        // Compute RY(A) = RY(B) + RY(C)
	        double RY_A = L.RY + R.RY;
	        
	        // Return new DPValues for the parent node A (parentDistance set to 0 since it's not specified)
	        return new DPValuesWeighted2(PX_A, PY_A, QX_A, QY_A, RX_A, RY_A, parentNode.getParentDistance());
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
	    private static final Map<BitSetPair, Double> weightCache = new ConcurrentHashMap<>();
		

		public CalculateWeightTask(STBipartition stb,
				ClusterCollection collection, List<Tree>trees) {
			this.stb = stb;
			this.containedClusterCollection = collection;
			this.trees = trees;
		}
		
		Double calculateMissingWeight(String method) {
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
		            return computeWeight(method);
		        });
	        } else {
	        	return computeWeight(method);
	        }
	        
		}
		
		Double getNc2(double n) {
			return (n*(n-1)*1.0)/2.0;
		}

		Double computeWeight(String method) {
//			System.out.println("\n\nhello\n\nMethod = " + method);
			//System.err.print("Calculating weight for: " + biggerSTB);
			Double weight = 0.0;
			
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
			
			
			
			Double w1 = 0.0;
			Double w2 = 0.0;
			
			// isNew variable controls whether the improved or the previous algorithm is being run
			// methods
			// 0 = unweighted
			// 1 = weighting with approach 1 - sum of three branch lengths
			// 2 - weighting with approach 2 - sum of two branch lengths
			
			boolean isTraversalImplementation = true;
//			isNew = true;
			
			// O(nk) traversal implementation instead of O(n^2k)
			if(isTraversalImplementation) {
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
		        
		        if(method.equals("base")) {
//		        	System.out.println("\n\n\nhi at base...\n\n\n");
					for (int t = 0; t < trees.size(); t++) {
						Tree tr = trees.get(t);
						
						Stack<DPValuesUnweighted> stack = new Stack<>();
						
						for (TNode node : tr.postTraverse()) {
							// System.err.println("Node is:" + node);
							if (node.isLeaf()) {
								String nodeName = node.getName();
								
								double PX, PY;
								
								
								// Check if the leaf matches with X or Y taxa
								if(xTaxaSet.contains(nodeName)) {
									PX = 1.0;
								} else {
									PX = 0.0;
								}
								
								if(yTaxaSet.contains(nodeName)) {
									PY = 1.0;
								} else {
									PY = 0.0;
								}
								
								DPValuesUnweighted mc = new DPValuesUnweighted(PX, PY);
				                
				                stack.push(mc);
							} else {
								cntt++;
								// For internal nodes, accumulate DPValues from children
								DPValuesUnweighted R = stack.pop();
								DPValuesUnweighted L = stack.pop();
			                	
								DPValuesUnweighted parentDPValues = getParentDPValues(L, R);
			                	
			                	stack.push(parentDPValues);
			                	
			                	double temp = 0.0;
			                	
	
			                	// Compute the weighted sum for each case
			                	double case1 = getNc2(L.PX) * R.PY;
			                	double case2 = L.PX * getNc2(R.PY);
			                	double case3 = getNc2(R.PX) * L.PY;
			                	double case4 = R.PX * getNc2(L.PY);
			                	
			                	// Update parentDPValues.QX and QY with the sum of the four cases
			                	weight += case1 + case2 + case3 + case4;
				        		
							}
						}
					}
				}
		        else if(method.equals("weighted_3_terminal")) {
//		        	System.out.println("\n\n\nat " + method + "...\n\n\n");
					for (int t = 0; t < trees.size(); t++) {
						Tree tr = trees.get(t);
						
						Stack<DPValuesWeighted1> stack = new Stack<>();
						
						for (TNode node : tr.postTraverse()) {
							// System.err.println("Node is:" + node);
							if (node.isLeaf()) {
								String nodeName = node.getName();
								
								double PX, PY, QX, QY;
								
								
								// Check if the leaf matches with X or Y taxa
								if(xTaxaSet.contains(nodeName)) {
									PX = 1.0;
								} else {
									PX = 0.0;
								}
								
								if(yTaxaSet.contains(nodeName)) {
									PY = 1.0;
								} else {
									PY = 0.0;
								}
								
								QX = 0.0;
								QY = 0.0;
								
								DPValuesWeighted1 mc = new DPValuesWeighted1(PX, PY, QX, QY, node.getParentDistance());
				                
				                stack.push(mc);
							} else {
								cntt++;
								// For internal nodes, accumulate DPValues from children
				                DPValuesWeighted1 R = stack.pop();
			                	DPValuesWeighted1 L = stack.pop();
			                	
			                	DPValuesWeighted1 parentDPValues = getParentDPValues(node, L, R);
			                	
			                	stack.push(parentDPValues);
			                	
			                	double temp = 0.0;
			                	
			                	// Compute the exponential decay factors
			                	double expWAB = Math.exp(-L.parentDistance); // e^(-w(A,B))
			                	double expWAC = Math.exp(-R.parentDistance); // e^(-w(A,C))
	
			                	// Compute the weighted sum for each case
			                	double case1 = expWAB * L.QX * (expWAC * R.PY); // e^(-w(A,B)) Q_X(B) e^(-w(A,C)) P_Y(C)
			                	double case2 = expWAB * L.PX * (expWAC * R.QY); // e^(-w(A,B)) P_X(B) e^(-w(A,C)) Q_Y(C)
			                	double case3 = expWAC * R.QX * (expWAB * L.PY); // e^(-w(A,C)) Q_X(C) e^(-w(A,B)) P_Y(B)
			                	double case4 = expWAC * R.PX * (expWAB * L.QY); // e^(-w(A,C)) P_X(C) e^(-w(A,B)) Q_Y(B)
			                	
			                	// Update parentDPValues.QX and QY with the sum of the four cases
			                	weight += case1 + case2 + case3 + case4;
				        		
							}
						}
					}
				}
				else if(method.equals("weighted_2_terminal")) {
					for (int t = 0; t < trees.size(); t++) {
						Tree tr = trees.get(t);
						
						Stack<DPValuesWeighted2> stack = new Stack<>();
						
						for (TNode node : tr.postTraverse()) {
							// System.err.println("Node is:" + node);
							if (node.isLeaf()) {
								String nodeName = node.getName();
								
								double PX, PY, QX, QY, RX, RY;
								
								
								// Check if the leaf matches with X or Y taxa
								if(xTaxaSet.contains(nodeName)) {
									PX = 1.0;
									RX = 1.0;
								} else {
									PX = 0.0;
									RX = 0.0;
								}
								
								if(yTaxaSet.contains(nodeName)) {
									PY = 1.0;
									RY = 1.0;
								} else {
									PY = 0.0;
									RY = 0.0;
								}
								
								QX = 0.0;
								QY = 0.0;
								
								DPValuesWeighted2 mc = new DPValuesWeighted2(PX, PY, QX, QY, RX, RY, node.getParentDistance());
				                
				                stack.push(mc);
							} else {
								cntt++;
								// For internal nodes, accumulate DPValues from children
				                DPValuesWeighted2 R = stack.pop();
			                	DPValuesWeighted2 L = stack.pop();
			                	
			                	DPValuesWeighted2 parentDPValues = getParentDPValues(node, L, R);
			                	
			                	stack.push(parentDPValues);
			                	
			                	double temp = 0.0;
			                	
			                	// Compute the weighted sum for each case
			                	double case1 = L.QX * R.RY; // Q_X(B) * R_Y(C)
			                	double case2 = L.RX * R.QY; // R_X(B) * Q_Y(C)
			                	double case3 = R.QX * L.RY; // Q_X(C) * R_Y(B)
			                	double case4 = R.RX * L.QY; // R_X(C) * Q_Y(B)
			                	
			                	// Update parentDPValues.QX and QY with the sum of the four cases
			                	weight += case1 + case2 + case3 + case4;
				        		
							}
						}
					}
				}
				
				w1 = weight;
				//System.out.println("Total "+cntt+" STBs");
			}
			else {   
				
				// the following is the previous implementation
				// Note that, this does one clever thing, rather than looping over all trees
				// it has a pre-determined aggregated list of all the unique STBs in the gene trees
				// basically it has a list of unique STBs and also a count map
				// so, it calculates for each of them separately, and multiplies by count of appearance of that STB, then adds to the total
				// since, practically there are many overlaps of STBs, this improves the performance a lot against the theoretical bound of O(n^2k)
				
				weight = 0.0;
				
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

		protected Double compute(String method) {
			return calculateMissingWeight(method);
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
