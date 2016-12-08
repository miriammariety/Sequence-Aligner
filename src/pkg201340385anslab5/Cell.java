/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pkg201340385anslab5;

/**
 *
 * @author ty
 */
public class Cell {
    Cell diagonal;
    Cell left;
    Cell top;
    int i;
    int j;
    int value;

    Cell(int value, int i, int j) {
        this.diagonal = null;
        this.left = null;
        this.top = null;
        this.i = i;
        this.j = j;
        this.value = value;
    }

    final int max(Cell diagonal, Cell left, Cell top, int score, int gap, int base) {
        int diagonalSum = diagonal.value + score;
        int leftSum = left.value + gap;
        int topSum = top.value + gap;
        int max = base;

        if (diagonalSum > max) {
            max = diagonalSum;
        }

        if (leftSum > max) {
            max = leftSum;
        }

        if (topSum > max) {
            max = topSum;
        }

        if (topSum == max) {
            this.top = top;
        }
        if (leftSum == max) {
            this.left = left;
        }
        if (diagonalSum == max) {
            this.diagonal = diagonal;
        }

        return max;
    }
}
