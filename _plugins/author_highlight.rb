# _plugins/author_highlight.rb
#
# Dynamically formats authors in bibliography citations using BibLaTeX annotation syntax:
#
# Highlighting:
#   AUTHOR+an = {1=highlight},
#   AUTHOR+an = {1=highlight; 2=highlight},
#
# Corresponding Author (adds superscript dagger †):
#   AUTHOR+an = {1=corresponding},
#   AUTHOR+an = {2=corresponding},
#   AUTHOR+an = {1=highlight, corresponding},
#   AUTHOR+an = {1=highlight; 2=corresponding},
#

module Jekyll
  module Filters
    def highlight_authors(reference, entry = nil)
      return reference unless reference && entry

      author_an = entry["author+an"] || entry[:"author+an"] || entry["author_an"]
      author_array = entry["author_array"]

      # Map of author index => { highlight: bool, corresponding: bool, author: Hash }
      annotated_authors = {}

      if author_an && !author_an.to_s.strip.empty?
        raw = author_an.to_s
        raw.scan(/(\d+)\s*[=:]\s*([^;]+)/).each do |idx_str, tag_str|
          idx = idx_str.to_i - 1
          tags = tag_str.downcase.split(/\s*,\s*/).map(&:strip)
          is_highlight = tags.any? { |t| t =~ /^(highlight|bold)$/i }
          is_corresponding = tags.any? { |t| t =~ /^(corresponding|corresp|corresponding-author)$/i }

          if (is_highlight || is_corresponding) && author_array && author_array[idx]
            annotated_authors[idx] = {
              highlight: is_highlight,
              corresponding: is_corresponding,
              author: author_array[idx]
            }
          end
        end
      end

      return reference if annotated_authors.empty?

      result = reference.dup

      annotated_authors.each do |_idx, info|
        name = info[:author]
        is_highlight = info[:highlight]
        is_corresponding = info[:corresponding]

        last = Regexp.escape(name["last"].to_s)
        first = name["first"].to_s
        initials = first.split(/[\s\-]+/).map { |part| "#{part[0]}." }.join(" ") rescue nil

        patterns = []
        if initials && !initials.empty?
          init_esc = Regexp.escape(initials)
          patterns << /(#{last},\s+#{init_esc})/
          patterns << /(#{init_esc}\s+#{last})/
        end

        if !first.empty?
          first_esc = Regexp.escape(first)
          patterns << /(#{last},\s+#{first_esc})/
          patterns << /(#{first_esc}\s+#{last})/
        end

        patterns << /(\b#{last},\s+[A-Z]\.)/

        patterns.each do |pat|
          if result =~ pat
            result = result.sub(pat) do |matched|
              formatted = matched
              formatted = "<b>#{formatted}</b>" if is_highlight
              formatted = "#{formatted}<sup>&dagger;</sup>" if is_corresponding
              formatted
            end
            break
          end
        end
      end

      result
    end
  end
end

# Ensure Liquid's Strainer registers the custom filter properly
if defined?(Liquid::Strainer)
  begin
    global_strainer = Liquid::Strainer.class_variable_get(:@@global_strainer)
    global_strainer.instance_variable_get(:@filter_methods).add("highlight_authors")
    Liquid::Strainer.class_variable_get(:@@strainer_class_cache).clear
    Liquid::Template.register_filter(Jekyll::Filters)
  rescue StandardError => _e
    Liquid::Template.register_filter(Jekyll::Filters)
  end
end
