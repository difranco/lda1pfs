package net.anthonydifranco.LDA;

import java.io.Serializable;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Map;
import java.util.List;

public class Encoder implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 4904891710882513758L;
	Map<Object, Integer> encode;
	List<Object> decode;
	int codesize;
	
	public Encoder() {
		encode = new TreeMap<Object, Integer>();
		decode = new ArrayList<Object>();
		codesize = 0;
	}

	public Integer encoded(Object arg) {
		if (! encode.containsKey(arg)) {
			encode.put(arg, codesize);
			decode.add(arg);
			codesize++;
		}
		return encode.get(arg);
	}

	public Object decoded(Integer arg) {
		return decode.get(arg);
	}
}
