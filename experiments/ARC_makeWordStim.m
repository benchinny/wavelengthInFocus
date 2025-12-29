% Download and install Optician Sans font to the OS:
% https://optician-sans.com/

%% set these parameters as needed

% Desired letter height in pixels
desiredHeight = 90;

% desired stroke width in pixels
desiredStroke = 4;

% text color (0-1 float)
color = [1 0 1];

% Number of 3-letter words we want to generate
numWords = 10;

% Letters we want to use
letters = {'a','c','e','m','n','o','p','q','r','s','u','v','w','x','y','z'};


%% sizing parameters - do not change

% compute desired ratio of stroke width to letter height, based on inputs
desiredRatio = desiredStroke/desiredHeight;

% Font
thisFont = 'Optician Sans';

% Image dimensions
imgHeight   = 1000;
imgWidth    = 2000;

% Initial font size
fontSize = 1000;

% Measured letter height in pixels at this font size (need to recalculate if font size changes)
measHeight = 500;

% Measured stroke width in pixels at this font size (need to recalculate if font size changes)
measStroke = 100;

% new stroke width needed to acheive desiredRatio at this font size
newStroke = desiredRatio*measHeight;

words = {'sea' 'sae' 'esa' 'eas' 'aes' 'ase'};

%% Generate 3 letter words
for x = 1:length(words)
    % Generate a random 3-letter word
    randLetters = letters(randperm(length(letters), 3));
    % word = strjoin(randLetters, '');
    word = words{x};

    % Create a figure (invisible)
    hFig = figure('Visible', 'off', 'Position', [100, 100, imgWidth, imgHeight], 'Color', 'k');
    hAx = axes('Parent', hFig, 'Position', [0, 0, 1, 1], 'Units', 'normalized', 'Color', 'k', 'XColor', 'none', 'YColor', 'none');

    % Set axis off
    axis off;

    % Create a text object
    hText = text(0.5, 0.5, word, 'FontName', thisFont, 'FontSize', fontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', color, 'Units', 'normalized');

    % Capture the frame and convert to image
    frame = getframe(hAx);
    img = frame2im(frame);

    % Close the figure
    close(hFig);

    % Convert the image to grayscale and then to binary
    grayImg = rgb2gray(img);
    binaryImg = imbinarize(grayImg);

    % Adjust stroke width using morphological operations
    if newStroke > measStroke
        % Increase stroke width using dilation
        se = strel('disk', round((newStroke - measStroke) / 2));
        adjustedImg = imdilate(binaryImg, se);
    else
        % Decrease stroke width using erosion
        se = strel('disk', round((measStroke - newStroke) / 2));
        adjustedImg = imerode(binaryImg, se);
    end

    % Convert binary image back to RGB
    finalImg = repmat(uint8(adjustedImg) * 255, 1, 1, 3);

    % Set the text color
    for c = 1:3
        channel = finalImg(:,:,c);
        channel(adjustedImg) = color(c) * 255;
        finalImg(:,:,c) = channel;
    end

    % Crop the black space around the letters
    stats = regionprops(adjustedImg, 'BoundingBox');
    if ~isempty(stats)
        bbox = stats(1).BoundingBox;
        for k = 2:numel(stats)
            bbox = [min(bbox(1), stats(k).BoundingBox(1)), ...
                min(bbox(2), stats(k).BoundingBox(2)), ...
                max(bbox(1) + bbox(3), stats(k).BoundingBox(1) + stats(k).BoundingBox(3)) - min(bbox(1), stats(k).BoundingBox(1)), ...
                max(bbox(2) + bbox(4), stats(k).BoundingBox(2) + stats(k).BoundingBox(4)) - min(bbox(2), stats(k).BoundingBox(2))];
        end
        croppedImg = imcrop(finalImg, bbox);
    else
        croppedImg = finalImg; % If no text is found, do not crop
    end

    % Resize to get letters at desired height
    croppedImg = imresize(croppedImg,desiredHeight/size(croppedImg,1));

    % Save the image
    folderName = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/';
    filename = sprintf([folderName 'word_image_%02d.png'], x);
    imwrite(croppedImg, filename);
end



