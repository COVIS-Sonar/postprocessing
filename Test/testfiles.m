

function [filename] = testfiles(mode,filetype)

  [filepath,name,ext]=fileparts(mfilename('fullpath'));

  root = fullfile(filepath,"../covis-test-data/old-covis-nas1/raw/2011/10/01");

  if filetype == "basename"
    switch mode
      case "imaging"
        filename = "APLUWCOVISMBSONAR001_20111001T210757.973Z-IMAGING";
      case "diffuse"
        filename = "APLUWCOVISMBSONAR001_20111001T215909.172Z-DIFFUSE";
      case "doppler"
        filename = "APLUWCOVISMBSONAR001_20130322T030922.697Z-DOPPLER";
    end
    return
  end

  filename = fullfile( root, testfiles(mode, "basename") );

  switch filetype
    case "gz"
      filename = strcat(filename,".tar.gz");
  end
