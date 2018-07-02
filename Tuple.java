package dspr_359_1pt688;

public class Tuple {

	public int x;
	public int y;
	
	public Tuple(){
		this.x = 0;
		this.y = 0;
	}
	
	public Tuple(int i, int k){
		this.x = i;
		this.y = k;
	}
	
	public void setX(int i){
		this.x = i;
	}
	public void setY(int i){
		this.y = i;
	}
	
	public int getX(){
		return this.x; 
	}
	
	public int getY(){
		return this.y;
	}
	
	public String toString(){
		return this.x + "," + this.y;
	}
	
	
	
}
