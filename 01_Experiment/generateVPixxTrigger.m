function VPixxRGBTrigger = generateVPixxTrigger(stimulusNumber, ditherProtection)
% generateVPixxTrigger Generates an RGB value that can be used as a trigger
% by displaying it as the top left pixel on a VPixx/EEG Monitor.
%
% VPixxRGBTrigger = generateVPixxTrigger(stimulusNumber, ditherProtection)
%
% --Inputs
%
% stimulusNumber, an interger whose value can be 0 to 255, inclusive.
% Any decimal value is rounded to an integer using MATLAB's round()
% function.
%
% ditherProtection, a boolean. TRUE: adds 8 to every color bit.
% Set to FALSE to disable this dither protection (You probably should set
% it to TRUE unless your specifically neeed to turn off dither protection 
% for some reason. If you don't know what dithering is, choose TRUE).
%
% -- Output
%
% VPixxRGBTrigger, a 1 x 3 matrix containing the R, G, and B values of the
% trigger.
%
% Coded: 02.2020, Suddha

assert(stimulusNumber >= 0 && stimulusNumber <= 255, 'Wrong input: stimulusNumber must be between 0 and 255 inclusive');
assert(islogical(ditherProtection), 'Wrong input: ditherProtection must be a logical (true/false)') 

stimulusNumber = uint8(round(stimulusNumber));

RByte = 0;
GByte = 0;
BByte = 0;

% RByte = bitand(stimulusNumber, 0b00001111) .* 2^4;
% GByte = bitand(stimulusNumber, 0b11110000);

RByte = bitand(stimulusNumber, bin2dec('00001111')) .* 2^4;
GByte = bitand(stimulusNumber, bin2dec('11110000'));

if(ditherProtection)
    RByte = RByte + 8;
    BByte = BByte + 8;
    GByte = GByte + 8;
end

VPixxRGBTrigger = [RByte GByte BByte];

