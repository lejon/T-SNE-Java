package com.jujutsu.tsne.barneshut;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.jujutsu.tsne.barneshut.VpTree.Node;

@SuppressWarnings("rawtypes")
public class TreePrinter {
	
	public interface AdditionalInfoProvider {
		public String provideInfo(Node node);
	}
	
	AdditionalInfoProvider provider;

	public TreePrinter() {
	}
	
	public TreePrinter(AdditionalInfoProvider additionalInfoProvider) {
		this.provider = additionalInfoProvider;
	}

	public <T extends Comparable<?>> void printNode(Node root) {
		int maxLevel = maxLevel(root);

		printNodeInternal(Collections.singletonList(root), 1, maxLevel);
	}

	private <T extends Comparable<?>> void printNodeInternal(List<Node> nodes, int level, int maxLevel) {
		if (nodes.isEmpty() || isAllElementsNull(nodes))
			return;

		int floor = maxLevel - level;
		int endgeLines = (int) Math.pow(2, (Math.max(floor - 1, 0)));
		int firstSpaces = (int) Math.pow(2, (floor)) - 1;
		int betweenSpaces = (int) Math.pow(2, (floor + 1)) - 1;

		printWhitespaces(firstSpaces);

		List<Node> newNodes = new ArrayList<Node>();
		for (Node node : nodes) {
			if (node != null) {
				System.out.print(node.index);
				if(provider!=null) {
					System.out.print("(" +provider.provideInfo(node) +")");
				}
				newNodes.add(node.getLeft());
				newNodes.add(node.getRight());
			} else {
				newNodes.add(null);
				newNodes.add(null);
				System.out.print(" ");
			}

			printWhitespaces(betweenSpaces);
		}
		System.out.println("");

		for (int i = 1; i <= endgeLines; i++) {
			for (int j = 0; j < nodes.size(); j++) {
				printWhitespaces(firstSpaces - i);
				if (nodes.get(j) == null) {
					printWhitespaces(endgeLines + endgeLines + i + 1);
					continue;
				}

				if (nodes.get(j).getLeft() != null)
					System.out.print("/");
				else
					printWhitespaces(1);

				printWhitespaces(i + i - 1);

				if (nodes.get(j).getRight() != null)
					System.out.print("\\");
				else
					printWhitespaces(1);

				printWhitespaces(endgeLines + endgeLines - i);
			}

			System.out.println("");
		}

		printNodeInternal(newNodes, level + 1, maxLevel);
	}

	private void printWhitespaces(int count) {
		for (int i = 0; i < count; i++)
			System.out.print(" ");
	}

	private <T extends Comparable<?>> int maxLevel(Node node) {
		if (node == null)
			return 0;

		return Math.max(maxLevel(node.getLeft()), maxLevel(node.getRight())) + 1;
	}

	private <T> boolean isAllElementsNull(List<T> list) {
		for (Object object : list) {
			if (object != null)
				return false;
		}

		return true;
	}


	public void printTreeHorizontal(Node node) {
		if (node.getRight() != null) {
			printTree(node.getRight(), true, "");
		}
		printNodeValue(node);
		if (node.getLeft() != null) {
			printTree(node.getLeft(), false, "");
		}
	}
	private void printNodeValue(Node node) {
		if (node == null) {
			System.out.print("<null>");
		} else {
			System.out.print(node.index);
			if(provider!=null) {
				System.out.print("(" +provider.provideInfo(node) +")");
			}
		}
		System.out.print('\n');
	}
	// use string and not stringbuffer on purpose as we need to change the indent at each recursion
	private void printTree(Node node, boolean isRight, String indent) {
		if (node.getRight() != null) {
			printTree(node.getRight(), true, indent + (isRight ? "        " : " |      "));
		}
		System.out.print(indent);
		if (isRight) {
			System.out.print(" /");
		} else {
			System.out.print(" \\");
		}
		System.out.print("----- ");
		printNodeValue(node);
		if (node.getLeft() != null) {
			printTree(node.getLeft(), false, indent + (isRight ? " |      " : "        "));
		}
	}

}
