function normalized_proj = normalize_projection(projection,bbox)

normalized_proj(1,:) = (projection(1,:) - bbox(1))/bbox(3);
normalized_proj(2,:) = (projection(2,:) - bbox(2))/bbox(4);
end