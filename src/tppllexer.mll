{
(** TODO The first special comment of the file is the comment associated with
    the whole module.*)
open Tpplparser

let keyword_table = Hashtbl.create 64
let _ =
  List.iter (fun (kwd, tok) -> Hashtbl.add keyword_table kwd tok)
    [

      (* Keywords *)
      "function",      FUNC;
      "if",            IF;
      "then",          THEN;
      "else",          ELSE;
      "utest",         UTEST;
      "observe",       OBSERVE;
      "match",         MATCH;
      "with",          WITH;

      (* Literals *)
      "true",          TRUE;
      "false",         FALSE;

      (* Symbolic Tokens *)
      "~",             TILDE;
      "(",             LPAREN;
      ")",             RPAREN;
      "{",             LCURLY;
      "}",             RCURLY;
      "[",             LSQUARE;
      "]",             RSQUARE;
      ":",             COLON;
      ";",             SEMICOLON;
      "::",            DCOLON;
      ",",             COMMA;
      ".",             DOT;
      "|",             VBAR;
      "->",            RARROW;

      (* Operators *)
      "=",             EQ;
      "+",             ADD;
      "-",             SUB;
      "*",             MUL;
      "/",             DIV;
      "%",             MOD;
      "<",             LESS;
      "<=",            LESSEQUAL;
      ">",             GREAT;
      ">=",            GREATEQUAL;
      "<<",            SHIFTLL;
      ">>",            SHIFTRL;
      ">>>",           SHIFTRA;
      "==",            EQUAL;
      "!=",            NOTEQUAL;
      "!",             NOT;
      "||",            OR;
      "&&",            AND;
      "++",            CONCAT;
    ]

}

let newline = '\n'
let whitespace = ['\t'' ']

let letter = ['A'-'Z''a'-'z']
let digit = ['0'-'9']
let nondigit = '_' | letter

let ident = (nondigit (digit | nondigit)*)

let integer = digit+
let real = digit* '.'? digit+ ( ['e''E'] ['+''-']? digit+ )?

let string_escape =
  "\\\\" | "\\\"" | "\\'" | "\\n" | "\\t" | "\\b" | "\\r" | "\\ "

let symtok =
  "~" | "(" | ")" | "{" | "}" | "[" | "]" | ":" | ";" | "::" |"," | "." | "|" |
  "->" | "=" | "+" | "-" | "*" | "/" | "%" | "<" | "<=" | ">" | ">=" | "<<" |
  ">>" | ">>>" | "==" | "!=" | "!" | "||" | "&&" | "++"

let line_comment = "//" [^'\n']*

rule main = parse
  | whitespace+ | line_comment { main lexbuf }
  | newline { Lexing.new_line lexbuf; main lexbuf }
  | "/*" { block_comment (Lexing.lexeme_start_p lexbuf) lexbuf; main lexbuf}
  | integer as str { INT(int_of_string str) }
  | real as str { FLOAT(float_of_string str) }
  | (ident as s) '('
      { match s with
        | "lam" -> LAM
        | _ -> FUNIDENT(s) }
  | ident | symtok as s
      { try Hashtbl.find keyword_table s with Not_found -> IDENT(s) }
  | '\'' (_ as c) '\'' { CHAR(c) }
  | '"' { read_string (Buffer.create 16) lexbuf }
  | eof { EOF }
  | _ { failwith "Unexpected char. TODO" }

and block_comment start_pos = parse
  | "*/" { }
  | eof { failwith "Unterminated block comment starting at TODO" }
  | newline { Lexing.new_line lexbuf; block_comment start_pos lexbuf }
  | _ { block_comment start_pos lexbuf }

and read_string buf = parse
  | '"' { STRING(Buffer.contents buf) }
  | string_escape as s
      { let c = match s with
          | "\\\\" -> '\\'
          | "\\\"" -> '"'
          | "\\'"  -> '\''
          | "\\n"  -> '\n'
          | "\\t"  -> '\t'
          | "\\b"  -> '\b'
          | "\\r"  -> '\r'
          | "\\ "  -> ' '
          | _ -> failwith "Should not happen"
        in Buffer.add_char buf c; read_string buf lexbuf }
  | '\\' newline whitespace*
    { Lexing.new_line lexbuf; read_string buf lexbuf }
  | '\\' { failwith "Invalid escape sequence (TODO)" }
  | eof { failwith "String not terminated (TODO)" }
  | [^ '"' '\\' '\n' ]+
      { Buffer.add_string buf (Lexing.lexeme lexbuf);
        read_string buf lexbuf }
  | newline { failwith "Unexpected line break in string. TODO" }
  | _ { failwith "Unexpected char in string. TODO" }

