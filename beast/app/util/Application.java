package beast.app.util;


import java.util.ArrayList;
import java.util.List;

import org.json.JSONObject;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;

@Description("BEAST application that handles argument parsing by introspection "
		+ "using Inputs declared in the class.")
public class Application {

	BEASTObject myBeastObject;

	public Application(BEASTObject myBeastObject) {
		this.myBeastObject = myBeastObject;
	}

	/** default input used for argument parsing **/
	protected Input<?> defaultInput = null;

	/**
	 * Arguments of the form -name value are processed by finding Inputs with
	 * matching name and setting their value.
	 * 
	 * If the input is a boolean that needs to be set to true, the 'value' a
	 * argument can be omitted.
	 * 
	 * The last argument is assigned to the defaultInput.
	 * **/
	public void parseArgs(String[] args, boolean sloppy) throws Exception {
		List<Input<?>> inputs = myBeastObject.listInputs();
		for (Input<?> input : inputs) {
			input.determineClass(this);
		}

		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			boolean done = false;
			if (arg.startsWith("-")) {
				String name = arg.substring(1);
				String value = (i < args.length - 1 ? args[i + 1] : null);
				for (Input<?> input : inputs) {
					if (input.getName().equals(name)) {
						try {
							if (input.getType() == Boolean.class) {
								if (value != null
										&& (value.toLowerCase().equals("true") || value
												.toLowerCase().equals("false"))) {
									input.setValue(value, null);
									i++;
								} else {
									input.setValue(Boolean.TRUE, null);
								}
							} else {
								input.setValue(value, null);
								i++;
							}
						} catch (Exception e) {
							throw new Exception("Problem parsing arguments:\n"
									+ e.getMessage());
						}
						done = true;
						break;
					}
				}
				if (name.equals("help")) {
					throw new Exception(""); // calling app should call getUsage()
				}
			} else {
				if (i == args.length - 1) {
					defaultInput.setValue(arg, null);
					done = true;
				}
			}
			if (!done) {
				if (sloppy) {
					Log.info.println("Unknown argument: " + args[i]
							+ " ignored.");
					i++;
				} else {
					throw new Exception("Unknown argument: " + args[i] + "\n");
							//+ getUsage());
				}
			}
		}

	}

	protected void parseArgs(JSONObject args) throws Exception {
		List<String> argList = new ArrayList<String>();
		for (String key : args.keySet()) {
			argList.add("-" + key.trim());
			argList.add(args.get(key).toString().trim());
		}
		parseArgs(argList.toArray(new String[] {}), true);
	}

	public String getUsage() {
		StringBuffer buf = new StringBuffer();
		try {
			List<Input<?>> inputs = myBeastObject.listInputs();
			buf.append("Usage: " + myBeastObject.getClass().getName() + "\n");
			for (Input<?> input : inputs) {
				buf.append("-" + input.getName() + " ");
				try {
					Class<?> typeClass = input.getType();
					if (typeClass == null) {
						input.determineClass(myBeastObject);
					}
					buf.append(input.getValueTipText());
				} catch (Exception e) {
					// ignore
				}
				buf.append("\t" + input.getTipText() + "\n");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		buf.append("-help\t show arguments");
		return buf.toString();
	}

	// placeholder, so that the main method compiles
	public void run() throws Exception {
		myBeastObject.initAndValidate();
	};

	// template for implementing a main for an application
	// the first argument is interpreted as class name of a BEASTObject
	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			BEASTObject myBeastObject = (BEASTObject) Class.forName(args[0])
					.newInstance();
			main = new Application(myBeastObject);
			String[] args2 = new String[args.length - 1];
			System.arraycopy(args, 1, args2, 0, args2.length);
			main.parseArgs(args2, false);
			main.run();
		} catch (Exception e) {
			System.out.println("Error:" + e.getMessage());
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
