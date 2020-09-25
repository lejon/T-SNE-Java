package com.jujutsu.tsne;
public class CheckResult {
		boolean check = true;
		String explanation = "";

		public CheckResult(boolean check, String explanation) {
			super();
			this.check = check;
			this.explanation = explanation;
		}
		
		public boolean check() {
			return check;
		}
		public void setCheck(boolean check) {
			this.check = check;
		}
		public String getExplanation() {
			return explanation;
		}
		public void setExplanation(String explanation) {
			this.explanation = explanation;
		} 
	}