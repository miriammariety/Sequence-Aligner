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
public class Sequence {

    // fasta sequence    
    protected String comment;
    protected String sequence;
    
    public Sequence (String comment, String sequence) {
        this.comment = comment;
        this.sequence = sequence.toUpperCase();
    }
    
    public boolean isProtein() {
        return sequence.matches("[ACDEFGHIKLMNPQRSTVWY]+") && 
                sequence.startsWith("M");
    }
    
    public boolean isNucleotide(){
        return sequence.matches("[actg]+");
    }
    
    @Override
    public String toString() {
        return this.comment + "\n" + this.sequence + "\n";
    }
    
    public String getID() {
        try {
            return this.comment.split("\\|")[1];
        } catch (Exception e) {
            return this.comment;
        }
    }
}

