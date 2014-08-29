package beast.app.tools;

import java.io.File;

import beast.app.beauti.Beauti;
import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;
import beast.app.util.ConsoleApp;
import beast.inference.PathSamplerFromFile;

//command line interface to PathSampler
public class PathSampler {
	
		
	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			PathSamplerFromFile sampler = new PathSamplerFromFile();
			
			if (args.length == 0) {
				// try the GUI version
				sampler.setID("PathSampler");
				sampler.initAndValidate();
				BeautiDoc doc = new BeautiDoc();
				doc.beautiConfig = new BeautiConfig();
				doc.beautiConfig.initAndValidate();
				doc.beautiConfig.suppressPlugins.add(sampler.getClass().getName() + ".mcmc");
				doc.beautiConfig.suppressPlugins.add(sampler.getClass().getName() + ".value");
				doc.beautiConfig.suppressPlugins.add(sampler.getClass().getName() + ".hosts");
				String fileSep = System.getProperty("file.separator");
				if (!sampler.model1Input.get().exists()) {
					sampler.model1Input.setValue(new File(Beauti.g_sDir + fileSep + "model.xml"), sampler);
				}
			
				BEASTObjectPanel panel = new BEASTObjectPanel(sampler, sampler.getClass(), doc);
				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
				if (dialog.showDialog()) {
					dialog.accept(sampler, doc);
					ConsoleApp app = new ConsoleApp("PathSampler", "Path Sampler: " + sampler.model1Input.get().getPath());
					sampler.initAndValidate();
					sampler.run();
				}
				return;
			}

			// continue with the command line version
			main = new Application(sampler);
			main.parseArgs(args, false);
			sampler.initAndValidate();
			sampler.run();
		} catch (Exception e) {
			//e.printStackTrace();
			System.out.println(e.getMessage());
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
