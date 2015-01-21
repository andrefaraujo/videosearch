import sys, string, os, datetime, time, glob, subprocess, re
from PIL import Image;

def parallelCommands(cmds, numThreads):
    if not cmds: return # empty list

    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)

    max_task = numThreads
    processes = []
    while True:
        while cmds and len(processes) < max_task:
            task = cmds.pop()
            print task
            processes.append(subprocess.Popen(task, shell=True))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            time.sleep(0.05)

currFolder = os.path.dirname(os.path.realpath(__file__));
rootFolder = "%s/../../../../videos/" % currFolder; #Change this path if necessary
overwriteExtractKeyframes = False;
overwriteExtractFeatures = False;
keyframeRate = 1;
keyframeHeight = 480;
numThreads = 5;
extractKeyframesScript = "python %s/../../indexer/extract_keyframes.py" % currFolder;
extractFeaturesScript = "%s/../../indexer/extract_features/extract" % currFolder;

months = ["201210", "201211", "201212", "201301"]; # light dataset
#months = ["201210", "201211", "201212", "201301", "201302", "201303", "201304", "201305", "201306", "201307", "201308", "201309"]; # full dataset
extractKeyframesCommands = [];
extractFeaturesCommands = [];

for month in months:
    dateFolders = glob.glob("%s/%s*" % (rootFolder, month));
    dateFolders.sort();
    for dateFolder in dateFolders:
        videoFiles = glob.glob("%s/*/*.700k.mp4" % dateFolder);
        videoFiles.sort();
        for videoFile in videoFiles:
            # Generate extract keyframes commands
            name, ext = os.path.splitext(videoFile);
            outputFolder = name + "_keyframes";
            if (not os.path.exists(outputFolder)) or overwriteExtractKeyframes:
                extractCmd = "%s %s %s %d scale=-1:%d" % (extractKeyframesScript, videoFile, outputFolder, keyframeRate, keyframeHeight);
                chmodCmd = "chmod -R a+rw %s" % outputFolder; #Change permissions of generated images
                extractKeyframesCommands.append(extractCmd + "; " + chmodCmd);
                continue;

            # Generate extract features commands
            keyframeFiles = glob.glob("%s/*.jpg" % outputFolder);
            featureFiles = glob.glob("%s/*.siftb" % outputFolder);
            if (len(keyframeFiles) != len(featureFiles)) or overwriteExtractFeatures:
                extractCmd = "%s -t %d -i %s.txt" % (extractFeaturesScript, numThreads, outputFolder);
                chmodCmd = "chmod -R a+rw %s" % outputFolder; #Change permissions of generated features
                extractFeaturesCommands.append(extractCmd + "; " + chmodCmd);
                continue;

    # Extract keyframes, convert keyframes
    parallelCommands(extractKeyframesCommands, numThreads);
    parallelCommands(extractFeaturesCommands, 1); # feature extraction already contains parallelization

