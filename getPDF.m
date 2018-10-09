function [output] = getPDF(y,mu,sigma)
%GETPDF 此处显示有关此函数的摘要
%   此处显示详细说明

    [M,N]=size(y);
    output=mvnpdf(y,mu,sigma);
end

